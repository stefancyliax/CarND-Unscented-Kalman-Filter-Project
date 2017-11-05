#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**********************************************
 * Initializes Unscented Kalman filter
 **********************************************/
UKF::UKF()
{
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);
  x_.fill(0.0);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  // A bicycle probably does not accelerate faster than 2m/s^2, so I'm choosing 1 m/s^2 here. (as shown in lesson 31)
  std_a_ = 1;

  // Process noise standard deviation yaw acceleration in rad/s^2
  // Lets say a bicycle turns at a rate of 60° per second at maximum. I'm choosing therefore PI/6. (as shown in lesson 31)
  std_yawdd_ = M_PI / 6;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**
  Complete the initialization. See ukf.h for other member properties.
  */

  // init flag to handle first measurement
  is_initialized_ = false;

  // init state dimension
  n_x_ = 5;
  // init augmented dimension
  n_aug_ = 7;

  // init previous timestep. Not really needed but for sanity reasons.
  time_us_ = 0;

  // define spreading parameter
  lambda_ = 3 - n_aug_;

  // use both sensor types
  use_laser_ = true;
  use_radar_ = true;

  // global parameter to control debugging output
  debugging_enabled_ = false;

  //set vector for weights
  weights_ = VectorXd(2 * n_aug_ + 1);
  double weight_0 = lambda_ / (lambda_ + n_aug_);
  weights_(0) = weight_0;
  for (int i = 1; i < 2 * n_aug_ + 1; i++)
  {
    double weight = 0.5 / (n_aug_ + lambda_);
    weights_(i) = weight;
  }

  // init lidar measurement noise covariance matrix
  R_lidar_ = MatrixXd(2, 2);
  R_lidar_ << std_laspx_ * std_laspx_, 0,
      0, std_laspy_ * std_laspy_;

  // init radar measurement noise covariance matrix
  R_radar_ = MatrixXd(3, 3);
  R_radar_ << std_radr_ * std_radr_, 0, 0,
      0, std_radphi_ * std_radphi_, 0,
      0, 0, std_radrd_ * std_radrd_;
}

// Destructor
UKF::~UKF() {}

/**********************************************
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 **********************************************/
void UKF::ProcessMeasurement(MeasurementPackage meas_package)
{

  // Initialize using the first measurement
  if (!is_initialized_)
  {
    UKF::InitializeState(meas_package);
    return;
  }

  // Calculate elapsed time since last measurement
  double dt = (meas_package.timestamp_ - time_us_) / 1000000.0; //dt - expressed in seconds
  if (debugging_enabled_)
    cout << "measurement received after " << dt << " seconds" << endl;
  time_us_ = meas_package.timestamp_;

  // Only do a prediction if the elapsed time is above 0.001s.
  // Otherwise just do another measurement update since the time is basically the same.
  if (dt > 0.001)
  {
    Prediction(dt);
  }

  // Distinguish between Radar and Lidar measurements
  // time of last update is only updated, if there really was a update
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_)
  {
    UpdateRadar(meas_package);
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_)
  {
    UpdateLidar(meas_package);
  }
}

/**********************************************
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 **********************************************/
void UKF::Prediction(double delta_t)
{

  //create sigma point matrix
  MatrixXd Xsig_aug_ = CreateSigmaPoints();

  // predict sigma points after delta_t
  Xsig_pred_ = PredictSigmaPoints(Xsig_aug_, delta_t);

  // predict state mean and state covariance matrix
  PredictMeanAndCovariance();
}

/**********************************************
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 **********************************************/
void UKF::UpdateLidar(MeasurementPackage meas_package)
{

  //set measurement dimension, lidar can measure px and py
  int n_z_ = 2;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig_ = MatrixXd(n_z_, 2 * n_aug_ + 1);

  //transform sigma points into measurement space
  //just using top 2 rows of the predicted sigma points,
  //instead of looping over the sigma points like in the Radar measurement
  Zsig_ = Xsig_pred_.topRows(2);

  //mean predicted measurement
  VectorXd z_pred_ = VectorXd(n_z_);
  z_pred_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  {
    z_pred_ = z_pred_ + weights_(i) * Zsig_.col(i);
  }

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z_, n_z_);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  { //2n+1 sigma points
    //residual
    VectorXd z_diff = Zsig_.col(i) - z_pred_;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  S = S + R_lidar_;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z_);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  { //2n+1 sigma points

    //residual
    VectorXd z_diff = Zsig_.col(i) - z_pred_;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3) > M_PI)
      x_diff(3) -= 2. * M_PI;
    while (x_diff(3) < -M_PI)
      x_diff(3) += 2. * M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  VectorXd z = VectorXd(n_z_);
  z = meas_package.raw_measurements_;
  //residual
  VectorXd z_diff = z - z_pred_;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();

  NIS_lidar_ = z_diff.transpose() * S.inverse() * z_diff;

  if (debugging_enabled_)
  {
    cout << "LIDAR input prediction " << endl
         << z_pred_ << endl;
    cout << "Sensor input " << endl
         << z << endl;
    cout << "Residual of LIDAR " << endl
         << z_diff << endl;
    cout << "update finished" << endl;
    cout << x_ << endl
         << endl;
  }
}

/**********************************************
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 **********************************************/
void UKF::UpdateRadar(MeasurementPackage meas_package)
{

  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z_ = 3;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig_ = MatrixXd(n_z_, 2 * n_aug_ + 1);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  { //2n+1 sigma points

    // extract values for better readability
    double p_x = Xsig_pred_(0, i);
    double p_y = Xsig_pred_(1, i);
    double v = Xsig_pred_(2, i);
    double yaw = Xsig_pred_(3, i);

    double v1 = cos(yaw) * v;
    double v2 = sin(yaw) * v;

    // measurement model
    float rho = sqrt(p_x * p_x + p_y * p_y);                             //rho
    float theta = atan2(p_y, p_x);                                       //theta
    float rho_dot = (p_x * v1 + p_y * v2) / sqrt(p_x * p_x + p_y * p_y); //rho_dot

    //handling if px and py are near zero
    if (p_x < 0.001 && p_y < 0.001)
    {
      theta = 0;
      rho_dot = 0;
    }

    Zsig_(0, i) = rho;
    Zsig_(1, i) = theta;
    Zsig_(2, i) = rho_dot;
  }

  //mean predicted measurement
  VectorXd z_pred_ = VectorXd(n_z_);
  z_pred_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  {
    z_pred_ = z_pred_ + weights_(i) * Zsig_.col(i);
  }

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z_, n_z_);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  { //2n+1 sigma points
    //residual
    VectorXd z_diff = Zsig_.col(i) - z_pred_;

    //angle normalization
    while (z_diff(1) > M_PI)
      z_diff(1) -= 2. * M_PI;
    while (z_diff(1) < -M_PI)
      z_diff(1) += 2. * M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix

  S = S + R_radar_;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z_);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  { //2n+1 sigma points

    //residual
    VectorXd z_diff = Zsig_.col(i) - z_pred_;
    //angle normalization
    while (z_diff(1) > M_PI)
      z_diff(1) -= 2. * M_PI;
    while (z_diff(1) < -M_PI)
      z_diff(1) += 2. * M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3) > M_PI)
      x_diff(3) -= 2. * M_PI;
    while (x_diff(3) < -M_PI)
      x_diff(3) += 2. * M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  VectorXd z = VectorXd(n_z_);
  z = meas_package.raw_measurements_;
  //residual
  VectorXd z_diff = z - z_pred_;

  //angle normalization
  while (z_diff(1) > M_PI)
    z_diff(1) -= 2. * M_PI;
  while (z_diff(1) < -M_PI)
    z_diff(1) += 2. * M_PI;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();

  NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;

  if (debugging_enabled_)
  {
    cout << "RADAR input prediction " << endl
         << z_pred_ << endl;
    cout << "Sensor input " << endl
         << z << endl;
    cout << "Residual of RADAR " << endl
         << z_diff << endl;
    cout << "update finished" << endl;
    cout << x_ << endl
         << endl;
  }
}

/**********************************************
 * Initialize the state when the first measurement arrives
 **********************************************/
//
void UKF::InitializeState(MeasurementPackage meas_package)
{
  double px = 0;
  double py = 0;

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
  {
    
    //Convert radar from polar to cartesian coordinates and initialize state.
    double rho = meas_package.raw_measurements_[0];
    double phi = meas_package.raw_measurements_[1];

    px = rho * cos(phi);
    py = rho * sin(phi);
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::LASER)
  {
    //Initialize state using LIDAR.
    px = meas_package.raw_measurements_[0];
    py = meas_package.raw_measurements_[1];
  }
  // initialize state vector
  // use first measurement for px and py.
  // in the sim the bicycle is moving, so we set the velocity to a sensible value for a bike (5m/s = 18 km/h).
  // from just the first measurement we can't really predict the direction the bike is moving. Therefore it is set to 0 (to the right) with a high uncertainty.
  x_ << px, py, 5, 0, 0;

  // initialize process covariance matrix
  // We are pretty certain about the position even after one measurement
  // After only one measurement we don't know anything about the velocity, so we are pretty uncertain.
  // We also don't know anything about the direction. Since we had to set it to a certain direction, we need a high uncertainty.
  P_ << 0.1, 0, 0, 0, 0,
      0, 0.1, 0, 0, 0,
      0, 0, 10, 0, 0,
      0, 0, 0, 10, 0,
      0, 0, 0, 0, 10;

  // save timestamp for use in next predict/measurement update
  time_us_ = meas_package.timestamp_;

  // print out initialisation
  if (debugging_enabled_)
    cout << "initialization finished" << endl;
  if (debugging_enabled_)
    cout << x_ << endl;

  // done initializing, no need to predict or update
  is_initialized_ = true;
}

/**********************************************
//Use current state to calculate augmented Sigma points
 **********************************************/
MatrixXd UKF::CreateSigmaPoints(void)
{
  //create augmented mean vector
  VectorXd x_aug_ = VectorXd(7);

  //create augmented state covariance
  MatrixXd P_aug_ = MatrixXd(7, 7);

  //create sigma point matrix
  MatrixXd Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  //create augmented mean state
  x_aug_.head(5) = x_;
  x_aug_(5) = 0;
  x_aug_(6) = 0;

  //create augmented covariance matrix
  P_aug_.fill(0.0);
  P_aug_.topLeftCorner(5, 5) = P_;
  P_aug_(5, 5) = std_a_ * std_a_;
  P_aug_(6, 6) = std_yawdd_ * std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug_.llt().matrixL();

  //create augmented sigma points
  Xsig_aug_.col(0) = x_aug_;
  for (int i = 0; i < n_aug_; i++)
  {
    Xsig_aug_.col(i + 1) = x_aug_ + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug_.col(i + 1 + n_aug_) = x_aug_ - sqrt(lambda_ + n_aug_) * L.col(i);
  }

  // print debug info
  if (debugging_enabled_)
  {
    cout << "Covariance augmented:" << endl
         << P_aug_ << endl;
    cout << "L Matrix: " << endl
         << L << endl;
    cout << "Sigma point : " << endl
         << Xsig_aug_ << endl;
  }

  return Xsig_aug_;
}

/**********************************************
 * Predict sigma points after delta_t 
 **********************************************/
MatrixXd UKF::PredictSigmaPoints(MatrixXd Xsig_aug_, double delta_t)
{
  //create matrix with predicted sigma points as columns
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  //predict sigma points
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  {
    //extract values for better readability
    double p_x = Xsig_aug_(0, i);
    double p_y = Xsig_aug_(1, i);
    double v = Xsig_aug_(2, i);
    double yaw = Xsig_aug_(3, i);
    double yawd = Xsig_aug_(4, i);
    double nu_a = Xsig_aug_(5, i);
    double nu_yawdd = Xsig_aug_(6, i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001)
    {
      px_p = p_x + v / yawd * (sin(yaw + yawd * delta_t) - sin(yaw));
      py_p = p_y + v / yawd * (cos(yaw) - cos(yaw + yawd * delta_t));
    }
    else
    {
      px_p = p_x + v * delta_t * cos(yaw);
      py_p = p_y + v * delta_t * sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd * delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5 * nu_a * delta_t * delta_t * cos(yaw);
    py_p = py_p + 0.5 * nu_a * delta_t * delta_t * sin(yaw);
    v_p = v_p + nu_a * delta_t;

    yaw_p = yaw_p + 0.5 * nu_yawdd * delta_t * delta_t;
    yawd_p = yawd_p + nu_yawdd * delta_t;

    //write predicted sigma point into right column
    Xsig_pred_(0, i) = px_p;
    Xsig_pred_(1, i) = py_p;
    Xsig_pred_(2, i) = v_p;
    Xsig_pred_(3, i) = yaw_p;
    Xsig_pred_(4, i) = yawd_p;
  }

  if (debugging_enabled_)
    cout << "Sigma point predicted: " << endl
         << Xsig_pred_ << endl;

  return Xsig_pred_;
}

/**********************************************
 * Predict state mean and state covariance
 **********************************************/
void UKF::PredictMeanAndCovariance(void)
{
  //predicted state mean
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  { //iterate over sigma points
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }

  //predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  { //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3) > M_PI)
      x_diff(3) -= 2. * M_PI;
    while (x_diff(3) < -M_PI)
      x_diff(3) += 2. * M_PI;

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
  }

  if (debugging_enabled_)
    cout << "Predicted state: " << endl
         << x_ << endl;
}
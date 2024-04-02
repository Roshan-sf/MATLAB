/* Roshan Jaiswal-Ferri Aero 220
Activity #9
Created February 8 2024
*/

// Setup pins and variables

// analog input pin for reading the potentiometer
// const int POT = A0;
// variable to capture pot value
// int potVal;

// Include servo library
#include <Servo.h>
#include <Adafruit_BNO055.h>

// PWM pin for controlling servo
const int SERVO = 9;

// value to command servo, in degrres; set to 90 degrees
int servoDeg = 90;

// instantiate servo
Servo myServo;
// instantiate sensor
Adafruit_BNO055 mySensor;

// create variables to load sensor data
imu::Vector<3> accel;
imu::Vector<3> gyro;
imu::Vector<3> mag;

void setup() {
  // set pin modes
  // pinMode(POT,INPUT);
  pinMode(SERVO,OUTPUT);
  // set initial pin states
  digitalWrite(SERVO,LOW);

  // start serial port, 9600 baud
  Serial.begin(9600);
  Serial.println("Serial started...");

  // Attach servo
  myServo.attach(SERVO);
  // set iniital position to 90 deg
  myServo.write(servoDeg);

  // start sensor
  mySensor.begin();
  // set sensor to AMG mode
  // mySensor.setMode(Adafruit_BNO055::OPERATION_MODE_AMG);
}

void loop() {
  // read potentiometer and report
  // potVal = analogRead(POT);
  // Serial.println("Pot value (ADC counts) is: ");
  // Serial.println(potVal);

  // read sensor data into variable, one vector per measurement type
  accel = mySensor.getVector(Adafruit_BNO055::VECTOR_ACCELEROMETER);
  gyro = mySensor.getVector(Adafruit_BNO055::VECTOR_GYROSCOPE);
  mag = mySensor.getVector(Adafruit_BNO055::VECTOR_MAGNETOMETER);

  // map pot value to servo position
  // servoDeg = map(potVal,0,1023,1,179);
  Serial.print("Servo position (degs) is: ");
  Serial.println(servoDeg);

  // command servo position based on mapped pot value
  myServo.write(servoDeg);

  // Output data to serial
  // Linear acceleration, m/sec^2
  Serial.print("\t\tAccel X: ");
  Serial.print(accel.x());
  Serial.print("\tY: ");
  Serial.print(accel.y());
  Serial.print("\tZ: ");
  Serial.print(accel.z());
  Serial.println("\tm/s^2");

  // move servo based on X axis acceleration
  servoDeg = map(accel.x(),-10,10,1,179);

  // set position
  myServo.write(servoDeg);

  delay(20);
} // end of loop()
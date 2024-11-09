#include <Wire.h>
#include <Servo.h>  // Include the built-in Arduino Servo library
 
Servo myServo;  // Create servo object to control a servo

// MPU-6050 I2C address
const int MPU_ADDR = 0x68;

// Global variables to store calibration offsets
float xOffset = 0.0;
float yOffset = 0.0;
float zOffset = 0.0;

// Proportional gain for the controller
float Kp = 1.5;  // You can adjust this value based on tuning

// Function to initialize the MPU-6050
void initializeMPU() {
    Wire.beginTransmission(MPU_ADDR);
    Wire.write(0x6B);  // Power management register
    Wire.write(0);     // Wake up the sensor (set sleep mode off)
    Wire.endTransmission(true);
    Serial.println("MPU-6050 initialization complete.");
}

// Function to calibrate the accelerometer
void calibrateAccelerometer(int numSamples) {
    long xSum = 0;
    long ySum = 0;
    long zSum = 0;

    // Gather specified number of samples for calibration
    for (int i = 0; i < numSamples; i++) {
        Wire.beginTransmission(MPU_ADDR);
        Wire.write(0x3B);
        Wire.endTransmission(false);
        Wire.requestFrom(MPU_ADDR, 6, true);

        if (Wire.available() == 6) {
            int16_t accX = Wire.read() << 8 | Wire.read();
            int16_t accY = Wire.read() << 8 | Wire.read();
            int16_t accZ = Wire.read() << 8 | Wire.read();

            xSum += accX;
            ySum += accY;
            zSum += accZ;
        }
        delay(10);
    }

    xOffset = xSum / numSamples;
    yOffset = ySum / numSamples;
    zOffset = (zSum / numSamples) - 16384;
    Serial.print("Calibration offsets - x: ");
    Serial.print(xOffset);
    Serial.print(", y: ");
    Serial.print(yOffset);
    Serial.print(", z: ");
    Serial.println(zOffset);
}

// Function to get roll angle
float getRollAngle() {
    Wire.beginTransmission(MPU_ADDR);
    Wire.write(0x3B);
    Wire.endTransmission(false);
    Wire.requestFrom(MPU_ADDR, 6, true);

    if (Wire.available() == 6) {
        int16_t rawX = Wire.read() << 8 | Wire.read();
        int16_t rawY = Wire.read() << 8 | Wire.read();
        int16_t rawZ = Wire.read() << 8 | Wire.read();

        float calibratedX = rawX - xOffset;
        float calibratedY = rawY - yOffset;
        float calibratedZ = rawZ - zOffset;

        float rollAngle = atan2(calibratedY, sqrt(calibratedX * calibratedX + calibratedZ * calibratedZ)) * (180.0 / PI);
        return rollAngle;
    } else {
        Serial.println("Error reading data");
        return 0;
    }
}

void setup() {
    Serial.begin(9600);
    Wire.begin();
    myServo.attach(9);  // Make sure pin 9 is correct
    initializeMPU();
    calibrateAccelerometer(100);

    // Testing initial servo response
    Serial.println("Testing servo movement...");
    myServo.write(0);   // Move to 0 degrees
    delay(1000);
    myServo.write(180); // Move to 180 degrees
    delay(1000);
    myServo.write(90);  // Move to 90 degrees (center)
    delay(1000);
    Serial.println("Servo test complete.");
}

void loop() {
    float rollAngle = getRollAngle();
    Serial.print("Roll Angle: ");
    Serial.println(rollAngle);

    float error = -rollAngle;
    float servoAngle = 90 + Kp * error;
    servoAngle = constrain(servoAngle, 0, 180);

    Serial.print("Servo Angle Command: ");
    Serial.println(servoAngle);

    myServo.write(servoAngle);
    delay(100);
}

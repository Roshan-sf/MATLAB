#include <Wire.h>

// Global variables to store calibration offsets
float xOffset = 0.0;
float yOffset = 0.0;
float zOffset = 0.0;

// Function to calibrate the accelerometer
void calibrateAccelerometer(int numSamples) {
    long xSum = 0;
    long ySum = 0;
    long zSum = 0;

    // Gather specified number of samples for calibration
    for (int i = 0; i < numSamples; i++) {
        Wire.beginTransmission(0x68);  // Assuming MPU-6050 address
        W ire.write(0x3B);              // Starting register for accelerometer readings
        Wire.endTransmission(false);
        Wire.requestFrom(0x68, 6, true);

        int16_t accX = Wire.read() << 8 | Wire.read();
        int16_t accY = Wire.read() << 8 | Wire.read();
        int16_t accZ = Wire.read() << 8 | Wire.read();

        xSum += accX;
        ySum += accY;
        zSum += accZ;

        delay(10);  // Small delay between samples
    }

    // Calculate the average values for offsets
    xOffset = xSum / numSamples;
    yOffset = ySum / numSamples;
    zOffset = (zSum / numSamples) - 16384;  // Assuming 1g (16384 raw value) on z-axis when flat

    // Print the calibration results to Serial Monitor
    Serial.print("Calibration offsets - x: ");
    Serial.print(xOffset);
    Serial.print(", y: ");
    Serial.print(yOffset);
    Serial.print(", z: ");
    Serial.println(zOffset);
}

// Function to get raw and calibrated readings from accelerometer
void getAccelerometerReadings() {
    // Start I2C transmission
    Wire.beginTransmission(0x68);  // MPU-6050 address
    Wire.write(0x3B);              // Starting register for accelerometer readings
    Wire.endTransmission(false);
    Wire.requestFrom(0x68, 6, true);

    // Read raw accelerometer data
    int16_t rawX = Wire.read() << 8 | Wire.read();
    int16_t rawY = Wire.read() << 8 | Wire.read();
    int16_t rawZ = Wire.read() << 8 | Wire.read();

    // Calibrate readings by subtracting offsets
    float calibratedX = rawX - xOffset;
    float calibratedY = rawY - yOffset;
    float calibratedZ = rawZ - zOffset;

    // Output the raw and calibrated readings to the Serial Monitor
    Serial.print("Raw - X: ");
    Serial.print(rawX);
    Serial.print(" Y: ");
    Serial.print(rawY);
    Serial.print(" Z: ");
    Serial.print(rawZ);

    Serial.print(" | Calibrated - X: ");
    Serial.print(calibratedX);
    Serial.print(" Y: ");
    Serial.print(calibratedY);
    Serial.print(" Z: ");
    Serial.println(calibratedZ);
}

void setup() {
    Serial.begin(9600);
    Wire.begin();

    // Calibrate the accelerometer using 100 samples
    calibrateAccelerometer(100);
}

void loop() {
    // Get and print accelerometer readings every 0.1 seconds (100 ms)
    getAccelerometerReadings();
    delay(100);
}

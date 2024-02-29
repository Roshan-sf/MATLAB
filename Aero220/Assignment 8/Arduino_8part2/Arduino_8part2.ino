/* Roshan Jaiswal-Ferri Aero 220
Activity #3 : PWM
Use PWM to vary LED brightness
Use potentiometer to vary PWM
Created February 8 2024
*/

// Motor pins; just "left side" of L293D
const int ENABLE1 = 11;
const int IN1 = 2;
const int IN2 = 4;

// Setup pins and variables

// analog input pin for reading the potentiometer
const int POT = A0;
// variable to capture pot value
int potVal = 0;

// Motor direction state; -1 -> reverse, 0 -> brake, 1 -> forward
int motorDir = 0;
// Motor PWM value "throttle"
int throttle = 0;

void setup() {
  // Setup pins
  pinMode(POT,INPUT);
  pinMode(ENABLE1,OUTPUT);
  // Force Motor 1 to DISABLE
  digitalWrite(ENABLE1,LOW);
  pinMode(IN1,OUTPUT);
  pinMode(IN2,OUTPUT);
}

void loop() {
  // Read Pot value
  potVal = analogRead(POT);
  //command motor based on potVal
  if (potVal < 502) {
    //reverse direction; check present dir first
    if (motorDir != -1) {
      // brake first
      brake();
      // set to reverse
      reverse();
    }
    throttle = map(potVal,502,0,0,255);
    speed(throttle);
  } else if (potVal > 522) {
    // forward direction; check present dir first
    if (motorDir != 1) {
      // brake first
      brake();
      // set to forward
      forward();
    }
    throttle = map(potVal,522,1023,0,255);
    speed(throttle);
  } else {
    // brake
    throttle = 0;
    brake();
  }
} // end of loop()

void brake() {
  // Stop rotation w/ brake
  digitalWrite(ENABLE1,LOW);
  digitalWrite(IN1,LOW);
  digitalWrite(IN2,LOW);
  digitalWrite(ENABLE1,HIGH);
  motorDir = 0;
}

void forward() {
  // Prepare to change Motor 1 state
  digitalWrite(ENABLE1,LOW);
  // Change state to "forward" Right-hand rotation
  digitalWrite(IN1,LOW);
  digitalWrite(IN2,HIGH);
  motorDir = 1;
}

void reverse() {
  // Prepare to change Motor 1 state
  digitalWrite(ENABLE1,LOW);
  // Change state to "reverse" Right-hand rotation
  digitalWrite(IN1,HIGH);
  digitalWrite(IN2,LOW);
  motorDir = -1;
}

void speed(int throttle) {
  analogWrite(ENABLE1,throttle);
}
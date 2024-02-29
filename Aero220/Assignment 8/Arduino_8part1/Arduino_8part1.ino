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

void setup() {
  // Setup pins
  pinMode(ENABLE1, OUTPUT);
  // Force Motor 1 to DISABLE
  digitalWrite(ENABLE1, LOW);
  pinMode(IN1, OUTPUT);
  pinMode(IN2, OUTPUT);
}

void loop() {
  // Start loop() with Motor 1 braked
  digitalWrite(ENABLE1,LOW);
  digitalWrite(IN1,LOW);
  digitalWrite(IN2,LOW);
  digitalWrite(ENABLE1,HIGH);
  delay(1000);

  // Prep to change Motor 1 state
  digitalWrite(ENABLE1,LOW);
  // Change state to "forward" Right-hand rotation
  digitalWrite(IN1,LOW);
  digitalWrite(IN2,HIGH);
  // Loop forward rotation by Enable PWM 0 -> 255
  for (int i = 0; i <= 255; i++) {
    // while spinning, leave IN1 & IN2 set and PWM ENABLE
    analogWrite(ENABLE1,i);
    delay(20);
  }

  // Stop rotation w/ brake
  digitalWrite(ENABLE1,LOW);
  digitalWrite(IN1,LOW);
  digitalWrite(IN2,LOW);
  digitalWrite(ENABLE1,HIGH);
  delay(1000);

  // Prep to change Motor 1 state
  digitalWrite(ENABLE1,LOW);
  // Change state to "reverse" Right-hand rotation
  digitalWrite(IN1,HIGH);
  digitalWrite(IN2,LOW);
  // Loop reverse rotation by Enable PWM 0 -> 255
  for (int i=0;i<=255;i++) {
    // while spinning, leave IN1 & IN2 set and PWM ENABLE
    analogWrite(ENABLE1,i);
    delay(20);
  }
} // end of loop()
/* Roshan Jaiswal-Ferri Aero 220
Activity #3 : PWM
Use PWM to vary LED brightness
Use potentiometer to vary PWM
Created February 8 2024
*/


const int LED = 11;
// potentiometer to A0 pin
const int POT = 0;
// variable to capture POT value
int potVal = 0;
// variable for the brightness
int bright = 0;

void setup() {
  // Configuring LED
  pinMode(LED,OUTPUT);
  digitalWrite(LED,LOW);
  Serial.begin(9600);
}

void loop() {
  // Vary LED brightness with POT:
  potVal = analogRead(POT);
  // map the input analog value to PWM
  bright = map((potVal), 0, 1023, 0, 255);
  analogWrite(LED,bright);

  Serial.print("Potentiometer Value: ");
  Serial.println(potVal);
  delay(1000);
  
}

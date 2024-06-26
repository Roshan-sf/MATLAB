/* Roshan Jaiswal-Ferri Aero 220
Activity #3 : PWM
Use PWM to vary LED brightness
Use potentiometer to vary PWM
Created February 8 2024
*/

const int LED = 11;


void setup() {
  // Configuring LED
  pinMode(LED,OUTPUT);
  digitalWrite(LED,LOW);
}

void loop() {
  // Vary LED brightness with PWM
  for (int i=0; i<=255; i++){
    analogWrite(LED,i);
    delay(10);
  }
// Make LED dimmer
  for (int i=255; i>=0; i--){
    analogWrite(LED,i);
    delay(10);
  }
}

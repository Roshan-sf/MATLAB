/* Roshan Jaiswal-Ferri Aero 220
Activity #1 : LED Blink
Will blink an external LED
Created January 25 2024
*/
const int LED = 13;
void setup() {
  // put your setup code here, to run once:
  pinMode(LED,OUTPUT);
}

void loop() {
  // put your main code here, to run repeatedly:
  for (int i=100; i<=1000; i+i+100)  {
    digitalWrite(LED,HIGH);
    delay(i);
    digitalWrite(LED,LOW);
    delay(i);
  }
}

/* Roshan Jaiswal-Ferri Aero 220
Activity #1 : Hello World Arduino
Will blink the on-board LED
Created January 11 2024
*/
//instead of using LED built in were assigning the LED variable to pin 13
const int LED = 13;

void setup() {
  // put your setup code here, to run once:
  //Setting the pin mode of the LED pin to an output
  pinMode(LED,OUTPUT);

}

void loop() {
  // put your main code here, to run repeatedly:
  for (int i=100; i<=1000; i=i+100)
  {
    //Turn LED on
    digitalWrite(LED,HIGH);
    //wait
    delay(i);
    //Turn LED off
    digitalWrite(LED,LOW);
    delay(i);
  }
}

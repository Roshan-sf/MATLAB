/* Roshan Jaiswal-Ferri Aero 220
Activity #4 : Input and Output
Will be able to Send & Recieve Digital & Analog Data
Created February 1 2024
*/

const int analogIn = 0;
float analogVal;
const int digitalOut = 12;
int digitalOutVal;
const int digitalIn = 11;
int digitalInVal;
const float adcScale = 5.0/1024.0;

void setup() {
  // Start serial on port, 9600:
  Serial.begin(9600);
  pinMode(analogIn,INPUT);
  pinMode(digitalOut,OUTPUT);
  pinMode(digitalIn,INPUT);

  digitalWrite(digitalOut, LOW);
  Serial.begin(9600);
}

void loop() {
  
  analogVal = analogRead(analogIn); // save analog reading
  digitalInVal = digitalRead(digitalIn); // save digital reading
    analogVal = analogVal * adcScale;
    if (analogVal < 2.5) {
      digitalWrite(digitalOut, HIGH);
    }
    else {
      digitalWrite(digitalOut, LOW);
    }
  Serial.print("Analog Input Value is: ");
  Serial.println(analogVal);
  Serial.print("Digital Input Value is: ");
  Serial.println(digitalInVal);
  delay(1000);
}

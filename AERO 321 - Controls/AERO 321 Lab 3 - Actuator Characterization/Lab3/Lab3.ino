
const int potPin = A3;   // Analog pin connected to potentiometer
const int servoPin = 2;  // Pin connected to the servo

void setup() {
  pinMode(servoPin, OUTPUT);
  Serial.begin(9600);
}

void loop() {
  // Read the analog value of the potentiometer
  int potValue = analogRead(potPin);
  // Map the raw ADC value to pulse length for the servo
  int pulseWidth = map(potValue, 0, 1023, 0, 3000); // 500-2500us typical range

  // Command the servo to move
  digitalWrite(servoPin, HIGH);
  delayMicroseconds(pulseWidth);
  digitalWrite(servoPin, LOW);
  delay(20 - pulseWidth / 1000); 

  // Print relevant infrormation
  Serial.print("ADC Value: ");
  Serial.print(potValue);
  Serial.print("  Pulse Width: ");
  Serial.println(pulseWidth);
}


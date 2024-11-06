// Set ultrasonic pins
const int trigPin = 9;
const int echoPin = 10;
// defines variables
const int motorPin1 = 2;
const int motorPin2 = 3;
const int enablePin = 5;

// distance function
float getDistance() {
  digitalWrite(trigPin, LOW);
  delayMicroseconds(2);
  digitalWrite(trigPin, HIGH);
  delayMicroseconds(10);
  digitalWrite(trigPin, LOW);

  long duration = pulseIn(echoPin, HIGH);
  float distance = (duration * 0.0343) / 2; // Convert to cm
  return distance;
}

//map motor speed
int mapDistanceToSpeed(float distance) {
  int speed = map((int)distance, 10, 100, 0, 255); // Map distance (10-100 cm) to speed (0-255)
  return constrain(speed, 0, 255); // Ensure speed is between 0 and 255
}

//motor pins
void setupMotor() {
  pinMode(motorPin1, OUTPUT);
  pinMode(motorPin2, OUTPUT);
  pinMode(enablePin, OUTPUT);
}


void move(int speed, int direction) {
  analogWrite(enablePin, speed);
  if (direction == 1) {
    digitalWrite(motorPin1, HIGH);
    digitalWrite(motorPin2, LOW);
  } else {
    digitalWrite(motorPin1, LOW);
    digitalWrite(motorPin2, HIGH);
  }
}

// stop function
void stopMotor() {
  digitalWrite(motorPin1, LOW);
  digitalWrite(motorPin2, LOW);
  analogWrite(enablePin, 0);
}

void setup() {
  Serial.begin(9600);
  setupMotor();
  pinMode(trigPin, OUTPUT);
  pinMode(echoPin, INPUT);
}

void loop() {
  float distance = getDistance();
  int speed = mapDistanceToSpeed(distance);
  move(speed, 1); // Move forward

  // Print distance and speed to Serial Monitor
  Serial.print(distance);
  Serial.print(" ");
  Serial.println(speed);

  delay(100); // Update every 0.1 seconds
}

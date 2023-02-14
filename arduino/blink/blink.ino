int LEDTime = 0;
static int LEDPin = 13;
int time1 = 0;
int time2 = 0;
int dur = 0;
float val = 0;
boolean started = false;


void setup() {
  // initialize digital pin LED_BUILTIN as an output.
    Serial.begin(57600);
  pinMode(LEDPin, OUTPUT);
    digitalWrite(LEDPin,LOW);
   time1 = millis();
}

// the loop function runs over and over again forever
void loop() {

   float val = Serial.parseFloat(); //read int or parseFloat for ..float...


   if (val == 1) { 
    started = true;
   }

   if (val == 0) { 
    started = false;
   }
   
   if started {
   time2 = millis();
   dur = time2-time1;

  if (dur > LEDTime) {
     digitalWrite(LEDPin,HIGH);
      delay(1000);
      digitalWrite(LEDPin,LOW);
      time1 = time2;
      
  LEDTime = random(2000, 3000);

  }

   Serial.println(dur);
   }
}

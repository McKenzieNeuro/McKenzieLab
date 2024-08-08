int LEDTime = 0;
static int LEDPin =13 ;
static int buttonPin = 8;
int time1 = 0;
int time2 = 0;
int time3 = 0;
int time4 = 0;
int dur = 0;
int dur_press = 0;


float val = 0;

boolean started = false;
boolean pressed = false;

void setup() {
  // initialize digital pin LED_BUILTIN as an output.
    Serial.begin(57600);
  pinMode(LEDPin, OUTPUT);
    digitalWrite(LEDPin,LOW);

     pinMode(buttonPin, INPUT);
      digitalWrite(buttonPin,HIGH);


    
   time1 = millis();
   time3 = millis();
}

// the loop function runs over and over again forever
void loop() {

   pressed = digitalRead(buttonPin); //read int or parseFloat for ..float...
 pressed = !pressed;
time4 = millis();
dur_press = time4 - time3;
 if (pressed == 1  & dur_press>5000 & !started) { 
 started = true;
  time3 = millis();
  dur_press = time4 - time3;

  Serial.println(started);
     

   
   }
   if (pressed == 1  & dur_press>5000 & started) { 
 started = false;
 time3 = millis();
  dur_press = time4 - time3;
   

   }
 
  
   if (started) {
   time2 = millis();
   dur = time2-time1;

  if (dur > LEDTime) {
     digitalWrite(LEDPin,HIGH);
      delay(200);
      digitalWrite(LEDPin,LOW);
      time1 = time2;
      
  LEDTime = random(2000, 3000);
   started = true;
  }

   //Serial.println(dur);
   }

   Serial.print(started);
   Serial.print(" ");
   Serial.print(pressed);
    Serial.print(" ");
   Serial.println(dur_press);
   
}

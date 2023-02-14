/*
  Blink

  Turns an LED on for one second, then off for one second, repeatedly.

  Most Arduinos have an on-board LED you can control. On the UNO, MEGA and ZERO
  it is attached to digital pin 13, on MKR1000 on pin 6. LED_BUILTIN is set to
  the correct LED pin independent of which board is used.
  If you want to know what pin the on-board LED is connected to on your Arduino
  model, check the Technical Specs of your board at:
  https://www.arduino.cc/en/Main/Products

  modified 8 May 2014
  by Scott Fitzgerald
  modified 2 Sep 2016
  by Arturo Guadalupi
  modified 8 Sep 2016
  by Colby Newman

  This example code is in the public domain.

  https://www.arduino.cc/en/Tutorial/BuiltInExamples/Blink
*/
int LEDTime = 0;
int sensorVal = 0;
boolean pushed = false;
int time1 = 0;
int time2 = 0;
int dur = 0;

// the setup function runs once when you press reset or power the board
void setup() {

  //start serial connection

  Serial.begin(9600);

  //configure pin 2 as an input and enable the internal pull-up resistor

  pinMode(2, INPUT_PULLUP);

  pinMode(13, OUTPUT);
   time1 = millis();

}

void loop() {

  //read the pushbutton value into a variable

  sensorVal = digitalRead(2);
 pushed = false;
  //print out the value of the pushbutton

  Serial.println(pushed);

  // Keep in mind the pull-up means the pushbutton's logic is inverted. It goes

  // HIGH when it's open, and LOW when it's pressed. Turn on pin 13 when the

  // button's pressed, and off when it's not:

  if (sensorVal == HIGH) {
  
    digitalWrite(13, LOW);

  } else {
    pushed = true;
  
  }
  
 
   if (pushed) {
    delay(500);
 sensorVal = digitalRead(2);
   }

   
  if (pushed) {
 
    while (sensorVal ==HIGH) {


  time2 = millis();
   dur = time2-time1;

  if (dur > LEDTime) {
     digitalWrite(13,HIGH);
      delay(1000);
      digitalWrite(13,LOW);
      time1 = time2;
      LEDTime = random(2000, 3000);
    
      sensorVal = digitalRead(2);



  

      
    }
    pushed = false;
     delay(2000);
  }

  
  }
}

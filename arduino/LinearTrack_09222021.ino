
static int solenoid1 = 13 ; // 
static int sensor1 = 9; // 


static int solenoid2 = 12 ; // 
static int sensor2 = 10; // 

int LEDout = 11;

static boolean ON = HIGH;
static boolean OFF = !ON;

int durL = 90; // duration left solenoid opens
int durR = 120; // duration right solenoid opens
int LED = 0;
int left = 0;     // variable for reading the pin status
int right = 0;
int time1 = 0;
int time2 = 0;
int dur = 0;

boolean middle = 0;
boolean donePulseM = false;
boolean donePulseL = false;
boolean donePulseR = false;

boolean nextLeft = 0;
boolean nextRight = 1;




void setup() {
      Serial.begin(9600);
      pinMode(solenoid1, OUTPUT);
      digitalWrite(solenoid1,LOW);

        pinMode(solenoid2, OUTPUT);
      digitalWrite(solenoid2,LOW);

       pinMode(LEDout, OUTPUT);
      digitalWrite(LEDout,LOW);
      
    
      pinMode(sensor1, INPUT);
     pinMode(sensor2, INPUT);
  digitalWrite(sensor1,HIGH);
   digitalWrite(sensor2,HIGH);
  time1 = millis();
  
}


void loop(){

  time2 = millis();
  dur = time2-time1;

  if (dur > LED) {
     digitalWrite(LEDout,HIGH);
      delay(100);
      digitalWrite(LEDout,LOW);
      time1 = time2;
      
  LED = random(20000, 30000);

  }

   left = digitalRead(sensor1);  // read input value
    right = digitalRead(sensor2);  // read input value

   Serial.print(right);

   Serial.print("  ");

   Serial.print(left);

    Serial.print("  ");

     Serial.println(dur);

   
    if (nextRight == true && !right ){
     digitalWrite(solenoid1,HIGH);
      delay(durR);
      digitalWrite(solenoid1,LOW);

  // Serial.println(right);
      nextRight = false;
      
      
    }

    if ( (nextRight ==false ) && (!left )) {
    digitalWrite(solenoid2,HIGH);
      delay(durL);
      digitalWrite(solenoid2,LOW);
 //    Serial.println(0);
 
      nextRight = true;
    }

}

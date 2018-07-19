/* 
 *  Immitation of Verasonics for Ultrafast Ultrasound Freehand Imaging.
 *  The Application is constantly waiting for triggering. Right After the trigger has happend, the system waits 2 seconds and send the trigger back.
 *  While Waiting Led on PIN 11 is activated
 *  
 */

#define LED     11
#define TRIGIN  9
#define TRIGOUT 13
#define DELAY   5*1e+3 // 5 sec
#define PWIDTH  200 //0.2 sec

static bool status;

void setup()
{
  pinMode(TRIGIN, INPUT);
  pinMode(TRIGOUT, OUTPUT);
  pinMode(LED, OUTPUT);

  digitalWrite(TRIGOUT, LOW);
  digitalWrite(LED, LOW);
  status = false;

}

void loop()
{
  status = digitalRead(TRIGIN);
  
  if(status){
    digitalWrite(LED, HIGH);
    delay(DELAY);
    digitalWrite(TRIGOUT, HIGH);
    digitalWrite(LED, LOW);
    delay(PWIDTH);
    digitalWrite(TRIGOUT, LOW);
  }

}

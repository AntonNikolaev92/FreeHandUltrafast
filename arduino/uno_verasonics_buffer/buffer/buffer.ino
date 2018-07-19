#define gpioVERIN 2
#define gpioBTNIN 4
#define gpioRPSTARTOUT 13
#define gpioRPSTOPOUT 12

#define PLENGTH 100

static bool interrupted;
static bool acquire;


void setup()
{

  pinMode(gpioVERIN, INPUT_PULLUP);
  pinMode(gpioBTNIN, INPUT);
  pinMode(gpioRPSTARTOUT, OUTPUT);
  pinMode(gpioRPSTOPOUT, OUTPUT);

  acquire = false;
  interrupted = false;

  attachInterrupt(digitalPinToInterrupt(gpioVERIN), procInterrupt, FALLING);
}

void loop()
{
  if (interrupted || digitalRead(gpioBTNIN)){
    if (acquire){
      digitalWrite(gpioRPSTARTOUT,HIGH);
      delay(PLENGTH);
      digitalWrite(gpioRPSTARTOUT,LOW);
    }else{
      digitalWrite(gpioRPSTOPOUT,HIGH);
      delay(PLENGTH);
      digitalWrite(gpioRPSTOPOUT,LOW);
    }
    acquire = !acquire;
    interrupted = false;
  }
}

void procInterrupt()
{
  interrupted = true;
}


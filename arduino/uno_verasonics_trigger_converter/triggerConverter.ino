static bool interrupted;

void setup()
{

  pinMode(2, INPUT_PULLUP);
  pinMode(13, OUTPUT);

  interrupted = false;

  attachInterrupt(digitalPinToInterrupt(2), procInterrupt, FALLING);
}

void loop()
{
  if (interrupted){
  digitalWrite(13,HIGH);
  delay(1000);
  digitalWrite(13,LOW);
  interrupted = false;
  }
}

void procInterrupt()
{
  interrupted = true;
}


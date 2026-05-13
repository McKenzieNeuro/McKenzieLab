/*
pulse_generator.ino

Arduino Uno pulse generator controlled by MATLAB

Command format:
pulseWidth_us,frequency_Hz,trainLength

Example:
200,1000,50
*/

const int pulsePin = 4;
const int LEDPin = 12;

String inputString = "";
bool commandComplete = false;

void setup()
{
    pinMode(pulsePin, OUTPUT);
    digitalWrite(pulsePin, LOW);

    Serial.begin(115200);

    inputString.reserve(50);

    Serial.println("READY");
}

void loop()
{
    if (commandComplete)
    {
        processCommand(inputString);

        inputString = "";
        commandComplete = false;
    }
}

void processCommand(String cmd)
{
    unsigned long pulseWidth;
    float frequency;
    int trainLength;

    int firstComma = cmd.indexOf(',');
    int secondComma = cmd.indexOf(',', firstComma + 1);

    if (firstComma < 0 || secondComma < 0)
    {
        Serial.println("ERR");
        return;
    }

    pulseWidth = cmd.substring(0, firstComma).toInt();

    frequency = cmd.substring(firstComma + 1, secondComma).toFloat();

    trainLength = cmd.substring(secondComma + 1).toInt();

    if (frequency <= 0)
    {
        Serial.println("ERR");
        return;
    }

    unsigned long period_ms =
        (unsigned long)(1000.0 / frequency);
    Serial.print("pulseWidth_ms: ");
    Serial.println(pulseWidth);
    Serial.print("Frequency: ");
    Serial.println(frequency);
    Serial.print("Train Length: ");
    Serial.println(trainLength);

    if (pulseWidth >= period_ms)
    {
        Serial.println("ERR");
        return;
    }

    generatePulseTrain(
        pulseWidth,
        period_ms,
        trainLength
    );

    Serial.println("DONE");
}

void generatePulseTrain(
    unsigned long pulseWidth,
    unsigned long period_ms,
    int trainLength)
{
    for (int i = 0; i < trainLength; i++)
    {
        digitalWrite(pulsePin, HIGH);

        delay(pulseWidth);

        digitalWrite(pulsePin, LOW);

        delay(period_ms - pulseWidth);
    }
}

void serialEvent()
{
    while (Serial.available())
    {
        char inChar = (char)Serial.read();

        if (inChar == '\n')
        {
            commandComplete = true;
        }
        else
        {
            inputString += inChar;
        }
    }
}

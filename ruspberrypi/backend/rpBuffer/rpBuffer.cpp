#include <stdio.h>
#include <wiringPi.h>
#include <cstring>
#include <stdlib.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <sys/select.h>
#include <pthread.h>
#include <limits.h>

#include "hidapi.h"

using namespace std;

/* Amfi Track HID*/
#define hidVID		0xc17
#define hidPID		0xd12

/* Ruspberry pi GPIO */
#define gpioIND		2 // 13
#define gpioACQSTOP		1 // 12
#define gpioACQSTART	0 // 11

#define szUSB		64 //bytes
#define nFRAMES_MAX		1000*szUSB
#define waitTIMEOUT     5000 // 5sec

#define PORT			5555

pthread_mutex_t mutex;
pthread_t idThread;

// Define buffer for coordinate reception 
unsigned char
	*bufRcvUSB, 		// usb data buffer
	*bufRcvPos;  // receive position data buffer
long *bufRcvTs;  // receive time stamp data buffer
unsigned int
	frameCur, 	// current displacement of writing pointer
	frameStart;
unsigned char *bufSend;

//networking connection variables
struct sockaddr_in server_addr;
int server;
fd_set waiting_set, waiting_set_init;

int errCode;

hid_device* hdev;

// timestamps
struct timespec ts1, ts0;

// support functions
int initialise(); 
int run();

int transferAll(int sock, unsigned char* buf, size_t szBuf);
void *acqDataThread( void* args );

int main(int argc, char *argv[])
{
	int status = 0;
	
	try
	{
		/* Initialization */
		status = initialise();
		digitalWrite(gpioIND, HIGH);
		delay(1000);
		digitalWrite(gpioIND, LOW);
		
		/* execution loop */
		while (!status)
			status = run();
	}
	catch (...)
	{

	}
	return 0;
}

/* Initialization section */
int initialise()
{
	mutex = PTHREAD_MUTEX_INITIALIZER;  // mutex for displacement variables	
		
	// RusberryPi GPIO
	wiringPiSetup();
	pinMode(gpioACQSTART, INPUT);
	pinMode(gpioACQSTOP, INPUT);
	pinMode(gpioIND, OUTPUT);
	digitalWrite(gpioIND, LOW);
	digitalWrite(gpioACQSTART, LOW);
	digitalWrite(gpioACQSTOP, LOW);

	// Memory
	bufRcvUSB = (unsigned char*)malloc(szUSB*sizeof(unsigned char));
	bufRcvPos = (unsigned char*)malloc(szUSB*nFRAMES_MAX*sizeof(unsigned char));
	bufSend = (unsigned char*)malloc( nFRAMES_MAX*(szUSB*sizeof(unsigned char) + sizeof(long)) );
	bufRcvTs = (long*)malloc(nFRAMES_MAX*sizeof(long));
	memset(bufRcvUSB, 0x00, szUSB*sizeof(unsigned char));
	memset(bufRcvPos, 0x00, szUSB*nFRAMES_MAX*sizeof(unsigned char));
	memset(bufRcvTs, 0x00, szUSB*sizeof(long));
	memset(bufSend, 0x00, nFRAMES_MAX*(szUSB*sizeof(unsigned char) + sizeof(long)));
	frameCur = 0;
	frameStart = 0;
	
	// AmfiTrack Device
	hdev = hid_open(hidVID, hidPID, NULL);
	if (!hdev) {
		printf("unable to open device\n");
		return 1;
	}
	hid_set_nonblocking(hdev, 0);
	
	int arg = 0;
	pthread_create(&idThread, NULL, acqDataThread, (void*)(&arg) );
	//pthread_join(idThread, NULL);
	
	// init server socket
	server = socket(AF_INET, SOCK_STREAM, 0);
	memset(&server_addr, '0', sizeof(server_addr)); 
	server_addr.sin_family = AF_INET;
	server_addr.sin_addr.s_addr = htonl(INADDR_ANY);
	server_addr.sin_port = htons(PORT);

	bind(server, (struct sockaddr*)&server_addr, sizeof(server_addr));
	
	FD_ZERO(&waiting_set);
	FD_SET(server, &waiting_set_init);
	listen(server, 1);
	// blind with indicator as everyng goes right
	
	digitalWrite(gpioIND, HIGH);
	delay(1000);
	digitalWrite(gpioIND, LOW);
	
	//pthread_join(idThread, NULL);
	return 0;
	
}

int run()
{
	// wait for the input triggering event
	while(!digitalRead(gpioACQSTART));
	pthread_mutex_lock(&mutex);
	frameStart = frameCur;
	pthread_mutex_unlock(&mutex);
	digitalWrite(gpioIND, HIGH);
	
	// wait for stop trigger input trigger
	int nFramesToSend = 0;
	unsigned char *ptrBufSend = bufSend;
	while (!digitalRead(gpioACQSTOP)) ;
	//usleep(1000);
	pthread_mutex_lock(&mutex);

	if (frameCur > frameStart){
		
		nFramesToSend = (frameCur - frameStart);
		memcpy(ptrBufSend, &nFramesToSend, sizeof(int));
		ptrBufSend += sizeof(int);
		memcpy( ptrBufSend, // positions
			bufRcvPos + frameStart*szUSB,
			(frameCur - frameStart)*szUSB*sizeof(unsigned char) );
		ptrBufSend += (frameCur - frameStart)*szUSB*sizeof(unsigned char);
		memcpy( ptrBufSend, // time stamps
			bufRcvTs  + frameStart,
			(frameCur - frameStart)*sizeof(long) );
	}
	else{
		nFramesToSend = (nFRAMES_MAX + frameCur - frameStart);
		memcpy( ptrBufSend, &nFramesToSend, sizeof(int) );
		ptrBufSend += sizeof(int);
		memcpy( ptrBufSend, // positions
			bufRcvPos + frameStart*szUSB,
			(nFRAMES_MAX - frameStart)*szUSB*sizeof(unsigned char) );
		ptrBufSend += (nFRAMES_MAX - frameStart)*szUSB*sizeof(unsigned char);
		memcpy( ptrBufSend,
			bufRcvPos,
			frameCur*szUSB*sizeof(unsigned char));
		ptrBufSend += frameCur*szUSB*sizeof(unsigned char);
		memcpy( ptrBufSend,// time stamp
			bufRcvTs  + frameStart,
			(nFRAMES_MAX - frameStart)*sizeof(long));
		ptrBufSend += (nFRAMES_MAX - frameStart)*sizeof(long);
		memcpy( ptrBufSend,
			bufRcvTs,
			frameCur*sizeof(long));
	}
		
	pthread_mutex_unlock(&mutex);
	digitalWrite(gpioIND, LOW);
	// remember last frame and swithc buffers
	
	// tcp connection and data transfer
	struct timeval timeout;
	timeout.tv_sec = 5;
	timeout.tv_usec = 0;
	waiting_set = waiting_set_init;
	int res = select(server + 1, &waiting_set, NULL, NULL, &timeout);
	if (res == 1) {
		int client = accept(server, NULL, NULL);
		transferAll( client, bufSend, nFramesToSend*(szUSB*sizeof(unsigned char) + sizeof(long) ) + sizeof(int) );
		printf(" Frames transfered: %i (%i bites)", nFramesToSend, nFramesToSend*(szUSB*sizeof(unsigned char) + sizeof(long)));
		close(client);
	}
	
	return 0;
}

/* keep sending data until all buffer is transfered */
int transferAll(int sock, unsigned char* buf, size_t szBuf)
{
	int numBytesTotal = 0;
	int numBytes = 0;
	int total = szBuf;
	
	while (numBytesTotal < szBuf){
		numBytes = write(sock, buf+numBytesTotal, szBuf - numBytesTotal);
		numBytesTotal += numBytes;
	}
	return 0;
}

/* Constantly acquire received data into the circular buffer*/
void *acqDataThread(void* args)
{
	
	struct timespec ts1, ts0;
	
	while (true) {
		hid_read(hdev, bufRcvUSB, szUSB);
		ts0 = ts1;
		clock_gettime(CLOCK_MONOTONIC, &ts1);
		pthread_mutex_lock(&mutex);
			if (ts1.tv_nsec > ts0.tv_nsec) bufRcvTs[frameCur] = (ts1.tv_nsec - ts0.tv_nsec);
			else bufRcvTs[frameCur] = (1e+9 + ts1.tv_nsec - ts0.tv_nsec);
			memcpy(bufRcvPos+frameCur*szUSB, bufRcvUSB, szUSB*sizeof(unsigned char));
			frameCur = frameCur < nFRAMES_MAX ? frameCur + 1 : 0;	
		pthread_mutex_unlock(&mutex);
		
	}
	
	return NULL;
}
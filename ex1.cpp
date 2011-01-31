#include <iostream>
#include <iomanip>
#include <windows.h>
#include <olectl.h>

#include <GL/glut.h>
#include <gl\gl.h>							// Header File For The OpenGL32 Library
#include <gl\glu.h>							// Header File For The GLu32 Library
#include <cv.h>
#include <highgui.h>

#include <math.h>

#include "PoseEstimation.h"
//#include "NeHeGL.h"
#pragma comment( lib, "opengl32.lib" )									// Search For OpenGL32.lib While Linking
#pragma comment( lib, "glu32.lib" )	

using namespace std;

// Added in Exercise 9 - Start *****************************************************************
#define PI 3.14159265

struct Position { double x,y,z; };

bool debugmode = false;
bool balldebug = false;
// Added in Exercise 9 - End *****************************************************************

int thresh = 100;
CvCapture* cap;

int bw_thresh = 100;

CvMemStorage* memStorage;

//left over from exercise 9
int towardscounter = 0;
Position ballpos;
int ballspeed = 100;

//matrices for the markers
float resultMatrix_005A[16]; // front cover
float resultMatrix_0272[16]; // image
float resultMatrix_1c44[16]; // turn right
float resultMatrix_1228[16]; // turn left

float resultTransposedMatrix[16];
float snowmanLookVector[4]; //use this for the model to marker calculation
int towards = 0x005A;
int towardsList[2] = {0x005a, 0x2720};

int markers[4] = {0x005a, 0x1228, 0x1c44, 0x0272};//marker database
GLuint	texture[5];	//number textures = number markers/number of all textures(for models,...)
//the bool to assign what markers were detected
bool detected[4];
int picture=3;
//1 = left, 2= right
int turn=0;
// Added in Exercise 9 - End *****************************************************************

//camera settings
const int width = 640; //important width, height must be correct!!!
const int height = 480;
const int camangle = 30;

unsigned char bkgnd[width*height*3];

void trackbarHandler(int pos) {
	thresh = pos;
}

void bw_trackbarHandler(int pos) {
	bw_thresh = pos;
}

//from NeHe, imports texture, binds it, returns the 
int BuildTexture( char* szPath, GLuint &texid)						// Load Image And Convert To A Texture
{
	HDC			hdcTemp;												// The DC To Hold Our Bitmap
	HBITMAP		hbmpTemp;												// Holds The Bitmap Temporarily
	IPicture	*pPicture;												// IPicture Interface
	OLECHAR		wszPath[MAX_PATH+1];									// Full Path To Picture (WCHAR)
	long		lWidth;													// Width In Logical Units
	long		lHeight;												// Height In Logical Units
	long		lWidthPixels;											// Width In Pixels
	long		lHeightPixels;											// Height In Pixels
	GLint		glMaxTexDim ;											// Holds Maximum Texture Size



	MultiByteToWideChar(CP_ACP, 0, szPath, -1, wszPath, MAX_PATH);		// Convert From ASCII To Unicode
	HRESULT hr = OleLoadPicturePath(wszPath, 0, 0, 0, IID_IPicture, (void**)&pPicture);

	if(FAILED(hr))														// If Loading Failed
		return FALSE;													// Return False

	hdcTemp = CreateCompatibleDC(GetDC(0));								// Create The Windows Compatible Device Context
	if(!hdcTemp)														// Did Creation Fail?
	{
		pPicture->Release();											// Decrements IPicture Reference Count
		return FALSE;													// Return False (Failure)
	}

	glGetIntegerv(GL_MAX_TEXTURE_SIZE, &glMaxTexDim);					// Get Maximum Texture Size Supported
	
	pPicture->get_Width(&lWidth);										// Get IPicture Width (Convert To Pixels)
	lWidthPixels	= MulDiv(lWidth, GetDeviceCaps(hdcTemp, LOGPIXELSX), 2540);
	pPicture->get_Height(&lHeight);										// Get IPicture Height (Convert To Pixels)
	lHeightPixels	= MulDiv(lHeight, GetDeviceCaps(hdcTemp, LOGPIXELSY), 2540);

	// Resize Image To Closest Power Of Two
	if (lWidthPixels <= glMaxTexDim) // Is Image Width Less Than Or Equal To Cards Limit
		lWidthPixels = 1 << (int)floor((log((double)lWidthPixels)/log(2.0f)) + 0.5f); 
	else  // Otherwise  Set Width To "Max Power Of Two" That The Card Can Handle
		lWidthPixels = glMaxTexDim;
 
	if (lHeightPixels <= glMaxTexDim) // Is Image Height Greater Than Cards Limit
		lHeightPixels = 1 << (int)floor((log((double)lHeightPixels)/log(2.0f)) + 0.5f);
	else  // Otherwise  Set Height To "Max Power Of Two" That The Card Can Handle
		lHeightPixels = glMaxTexDim;
	
	//	Create A Temporary Bitmap
	BITMAPINFO	bi = {0};												// The Type Of Bitmap We Request
	DWORD		*pBits = 0;												// Pointer To The Bitmap Bits

	bi.bmiHeader.biSize			= sizeof(BITMAPINFOHEADER);				// Set Structure Size
	bi.bmiHeader.biBitCount		= 32;									// 32 Bit
	bi.bmiHeader.biWidth		= lWidthPixels;							// Power Of Two Width
	bi.bmiHeader.biHeight		= lHeightPixels;						// Make Image Top Up (Positive Y-Axis)
	bi.bmiHeader.biCompression	= BI_RGB;								// RGB Encoding
	bi.bmiHeader.biPlanes		= 1;									// 1 Bitplane

	//	Creating A Bitmap This Way Allows Us To Specify Color Depth And Gives Us Imediate Access To The Bits
	hbmpTemp = CreateDIBSection(hdcTemp, &bi, DIB_RGB_COLORS, (void**)&pBits, 0, 0);
	
	if(!hbmpTemp)														// Did Creation Fail?
	{
		DeleteDC(hdcTemp);												// Delete The Device Context
		pPicture->Release();											// Decrements IPicture Reference Count
		return FALSE;													// Return False (Failure)
	}

	SelectObject(hdcTemp, hbmpTemp);									// Select Handle To Our Temp DC And Our Temp Bitmap Object

	// Render The IPicture On To The Bitmap
	pPicture->Render(hdcTemp, 0, 0, lWidthPixels, lHeightPixels, 0, lHeight, lWidth, -lHeight, 0);

	// Convert From BGR To RGB Format And Add An Alpha Value Of 255
	for(long i = 0; i < lWidthPixels * lHeightPixels; i++)				// Loop Through All Of The Pixels
	{
		BYTE* pPixel	= (BYTE*)(&pBits[i]);							// Grab The Current Pixel
		BYTE  temp		= pPixel[0];									// Store 1st Color In Temp Variable (Blue)
		pPixel[0]		= pPixel[2];									// Move Red Value To Correct Position (1st)
		pPixel[2]		= temp;											// Move Temp Value To Correct Blue Position (3rd)

		// This Will Make Any Black Pixels, Completely Transparent		(You Can Hardcode The Value If You Wish)
		if ((pPixel[0]==0) && (pPixel[1]==0) && (pPixel[2]==0))			// Is Pixel Completely Black
			pPixel[3]	=   0;											// Set The Alpha Value To 0
		else															// Otherwise
			pPixel[3]	= 255;											// Set The Alpha Value To 255
	}

	glGenTextures(1, &texid);											// Create The Texture

	// Typical Texture Generation Using Data From The Bitmap
	glBindTexture(GL_TEXTURE_2D, texid);								// Bind To The Texture ID
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);		// (Modify This For The Type Of Filtering You Want)
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);     // (Modify This For The Type Of Filtering You Want)
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, lWidthPixels, lHeightPixels, 0, GL_RGBA, GL_UNSIGNED_BYTE, pBits);	// (Modify This If You Want Mipmaps)

	DeleteObject(hbmpTemp);												// Delete The Object
	DeleteDC(hdcTemp);													// Delete The Device Context

	pPicture->Release();												// Decrements IPicture Reference Count

	return TRUE;														// Return True (All Good)
}

void initVideoStream() {
	cap = cvCaptureFromCAM (0);

	if (!cap) {
		cout << "No webcam found, using video file\n";
		cap = cvCaptureFromFile("C:\\Andi\\Uni\\SVNs\\far_intern\\teaching\\2010WS\\AR\\Exercises\\Solutions\\MarkerMovie.mpg");
		if (!cap) {
			cout << "No video file found. Exiting.\n";
			exit(0);
		}
	}
	//bind textures here, ALL TEXTURES!
	BuildTexture("C:\\Users\\MRX\\Documents\\Visual Studio 2010\\Projects\\ex1\\Debug\\funny-pictures-cat-goes-pew.jpg", texture[0]);
	BuildTexture("C:\\Users\\MRX\\Documents\\Visual Studio 2010\\Projects\\ex1\\chinese_curry.jpg", texture[1]);
	BuildTexture("C:\\Users\\MRX\\Documents\\Visual Studio 2010\\Projects\\ex1\\toybox01_089.jpg", texture[2]);
	BuildTexture("C:\\Users\\MRX\\Documents\\Visual Studio 2010\\Projects\\ex1\\toybox01_092.jpg", texture[3]);
	BuildTexture("C:\\Users\\MRX\\Documents\\Visual Studio 2010\\Projects\\ex1\\Debug\\funny-pictures-cat-goes-pew.jpg", texture[4]);
}

int subpixSampleSafe ( const IplImage* pSrc, CvPoint2D32f p )
{
	int x = int( floorf ( p.x ) );
	int y = int( floorf ( p.y ) );

	if ( x < 0 || x >= pSrc->width  - 1 ||
		 y < 0 || y >= pSrc->height - 1 )
		return 127;

	int dx = int ( 256 * ( p.x - floorf ( p.x ) ) );
	int dy = int ( 256 * ( p.y - floorf ( p.y ) ) );

	unsigned char* i = ( unsigned char* ) ( ( pSrc->imageData + y * pSrc->widthStep ) + x );
	int a = i[ 0 ] + ( ( dx * ( i[ 1 ] - i[ 0 ] ) ) >> 8 );
	i += pSrc->widthStep;
	int b = i[ 0 ] + ( ( dx * ( i[ 1 ] - i[ 0 ] ) ) >> 8 );
	return a + ( ( dy * ( b - a) ) >> 8 );
}

// Added in Exercise 9 - Start *****************************************************************
void multMatrix(float mat[16], float vec[4])
{
	for(int i=0; i<4; i++)
	{
		snowmanLookVector[i] = 0;
		for(int j=0; j<4; j++)
			  snowmanLookVector[i] += mat[4*i + j] * vec[j];
	}
}

void moveBall(float mat[16])
{
	float vector[3];
	vector[0] = mat[3] - ballpos.x;
	vector[1] = mat[7] - ballpos.y;
	vector[2] = mat[11] - ballpos.z;

	float length = sqrt( vector[0]*vector[0] + vector[1]*vector[1] + vector[2]*vector[2] );
	if(balldebug) std::cout << length << std::endl;
	if( length < 0.01) 
	{ 
		towards = towardsList[(towardscounter++)%2]; 
		if(balldebug) std::cout << "target changed to marker " << towards << std::endl; 
		ballspeed = 60 + 80 * rand()/RAND_MAX;
		return; 
	}
	ballpos.x += vector[0] / (ballspeed * length);
	ballpos.y += vector[1] / (ballspeed * length);
	ballpos.z += vector[2] / (ballspeed * length);

}

void rotateToMarker(float thisMarker[16], float lookAtMarker[16], int markernumber)
{
	float vector[3];
	vector[0] = lookAtMarker[3] - thisMarker[3];
	vector[1] = lookAtMarker[7] - thisMarker[7];
	vector[2] = lookAtMarker[11] - thisMarker[11];

	if(towards == markernumber) moveBall(lookAtMarker);
	
	//normalize vector
	float help = sqrt( vector[0]*vector[0] + vector[1]*vector[1] + vector[2]*vector[2] );
	vector[0] /= help;
	vector[1] /= help;
	vector[2] /= help;

	if(debugmode) std::cout << "Vector: " << vector[0] << ", " << vector[1] << ", " << vector[2] << std::endl;

	float defaultLook[4] = {1,0,0,0};
	multMatrix(thisMarker, defaultLook);

	//normalize snowmanLookVector
	help = sqrt( snowmanLookVector[0]*snowmanLookVector[0] + snowmanLookVector[1]*snowmanLookVector[1] + snowmanLookVector[2]*snowmanLookVector[2] );
	snowmanLookVector[0] /= help;
	snowmanLookVector[1] /= help;
	snowmanLookVector[2] /= help;

	if(debugmode) std::cout << "SnowmanLookVector: " << snowmanLookVector[0] << ", " << snowmanLookVector[1] << ", " << snowmanLookVector[2] << std::endl;

	float angle = (180 / PI) * acos( vector[0] * snowmanLookVector[0] + vector[1] * snowmanLookVector[1] + vector[2] * snowmanLookVector[2] );
	if((vector[0] * snowmanLookVector[1] - vector[1] * snowmanLookVector[0]) < 0 ) angle *= -1;
	
	if(debugmode) std::cout << "Angle: " << angle << std::endl;
	
	glRotatef(angle, 0, 0, 1);
}
// Added in Exercise 9 - End *****************************************************************

void init()
{
	cvNamedWindow ("Exercise 8 - Original Image", CV_WINDOW_AUTOSIZE);
	cvNamedWindow ("Exercise 8 - Converted Image", CV_WINDOW_AUTOSIZE);
	cvNamedWindow ("Exercise 8 - Stripe", CV_WINDOW_AUTOSIZE);
	cvNamedWindow ("Marker", 0 );
	cvResizeWindow("Marker", 120, 120 );
	initVideoStream();

	int value = thresh;
	int max = 255;
	cvCreateTrackbar( "Threshold", "Exercise 8 - Converted Image", &value, max, trackbarHandler);

	int bw_value = bw_thresh;
	cvCreateTrackbar( "BW Threshold", "Exercise 8 - Converted Image", &bw_value, max, bw_trackbarHandler);

	memStorage = cvCreateMemStorage();
}
//calculate deteted markers
void idle()
{
	detected[0]=detected[1]=detected[2]=detected[3]=false;
	bool isFirstStripe = true;

	bool isFirstMarker = true;

	IplImage* iplGrabbed = cvQueryFrame(cap);

	if(!iplGrabbed){
		printf("Could not query frame. Trying to reinitialize.\n");
		cvReleaseCapture (&cap);
		initVideoStream();
		return;
	}

	CvSize picSize = cvGetSize(iplGrabbed);

	memcpy( bkgnd, iplGrabbed->imageData, sizeof(bkgnd) );

	IplImage* iplConverted = cvCreateImage(picSize, IPL_DEPTH_8U, 1);
	IplImage* iplThreshold = cvCreateImage(picSize, IPL_DEPTH_8U, 1);

	cvConvertImage(iplGrabbed, iplConverted, 0);
	cvThreshold(iplConverted, iplThreshold, thresh, 255, CV_THRESH_BINARY);
	//cvAdaptiveThreshold(iplConverted, iplThreshold, 255, CV_ADAPTIVE_THRESH_MEAN_C, CV_THRESH_BINARY, 33, 5);

	// Find Contours
	CvSeq* contours;

	cvFindContours(
		iplThreshold, memStorage, &contours, sizeof(CvContour),
		CV_RETR_LIST, CV_CHAIN_APPROX_SIMPLE
	);

	for (; contours; contours = contours->h_next)
	{
		CvSeq* result = cvApproxPoly(
			contours, sizeof(CvContour), memStorage, CV_POLY_APPROX_DP,
			cvContourPerimeter(contours)*0.02, 0
		);

		CvRect r = cvBoundingRect(result);
		if (r.height < 20 || r.width < 20 || r.height >= iplGrabbed->height - 10 || r.width >= iplGrabbed->width - 10) {
			continue;
		}

		if (result->total==4)
		{
			int count = 4;
			CvPoint *rect = new CvPoint[4];
			cvCvtSeqToArray(result, rect);
			cvPolyLine(iplGrabbed, &rect, &count, 1, 1, CV_RGB(255,0,0), 2);
			
			float lineParams[16];

			for (int i=0; i<4; ++i)
			{
				cvCircle (iplGrabbed, rect[i], 3, CV_RGB(0,255,0), -1);

				double dx = (double)(rect[(i+1)%4].x-rect[i].x)/7.0;
				double dy = (double)(rect[(i+1)%4].y-rect[i].y)/7.0;

				int stripeLength = (int)(0.8*sqrt (dx*dx+dy*dy));
				if (stripeLength < 5)
				stripeLength = 5;

				//make stripeLength odd (because of the shift in nStop)
				stripeLength |= 1;

				//e.g. stripeLength = 5 --> from -2 to 2
				int nStop  = stripeLength>>1;
				int nStart = -nStop;

				CvSize stripeSize;
				stripeSize.width = 3;
				stripeSize.height = stripeLength;

				CvPoint2D32f stripeVecX;
				CvPoint2D32f stripeVecY;

				//normalize vectors
				double diffLength = sqrt ( dx*dx+dy*dy );
				stripeVecX.x = dx / diffLength;
				stripeVecX.y = dy / diffLength;

				stripeVecY.x =  stripeVecX.y;
				stripeVecY.y = -stripeVecX.x;

				IplImage* iplStripe = cvCreateImage( stripeSize, IPL_DEPTH_8U, 1 );

				// Array for edge point centers
				CvPoint2D32f points[6];

				for (int j=1; j<7; ++j)
				{
					double px = (double)rect[i].x+(double)j*dx;
					double py = (double)rect[i].y+(double)j*dy;

					CvPoint p;
					p.x = (int)px;
					p.y = (int)py;
					cvCircle (iplGrabbed, p, 2, CV_RGB(0,0,255), -1);

					for ( int m = -1; m <= 1; ++m )
					{
						for ( int n = nStart; n <= nStop; ++n )
						{
							CvPoint2D32f subPixel;

							subPixel.x = (double)p.x + ((double)m * stripeVecX.x) + ((double)n * stripeVecY.x);
							subPixel.y = (double)p.y + ((double)m * stripeVecX.y) + ((double)n * stripeVecY.y);

							CvPoint p2;
							p2.x = (int)subPixel.x;
							p2.y = (int)subPixel.y;

							if (isFirstStripe)
								cvCircle (iplGrabbed, p2, 1, CV_RGB(255,0,255), -1);
							else
								cvCircle (iplGrabbed, p2, 1, CV_RGB(0,255,255), -1);

							int pixel = subpixSampleSafe (iplConverted, subPixel);

							int w = m + 1; //add 1 to shift to 0..2
							int h = n + ( stripeLength >> 1 ); //add stripelenght>>1 to shift to 0..stripeLength

							*(iplStripe->imageData + h * iplStripe->widthStep  + w) =  pixel; //set pointer to correct position and safe subpixel intensity
						}
					}

					//use sobel operator on stripe
					// ( -1 , -2, -1 )
					// (  0 ,  0,  0 )
					// (  1 ,  2,  1 )
					
					double* sobelValues = new double[stripeLength-2];
					for (int n = 1; n < (stripeLength-1); n++)
					{
						unsigned char* stripePtr = ( unsigned char* )( iplStripe->imageData + (n-1) * iplStripe->widthStep );
						double r1 = -stripePtr[ 0 ] - 2 * stripePtr[ 1 ] - stripePtr[ 2 ];

						stripePtr += 2*iplStripe->widthStep;
						double r3 =  stripePtr[ 0 ] + 2 * stripePtr[ 1 ] + stripePtr[ 2 ];
						sobelValues[n-1] = r1+r3;
					}

					double maxVal = -1;
					int maxIndex = 0;
					for (int n=0; n<stripeLength-2; ++n)
					{
						if ( sobelValues[n] > maxVal )
						{
							maxVal = sobelValues[n];
							maxIndex = n;
						}
					}

					double y0,y1,y2; // y0 .. y1 .. y2
					y0 = (maxIndex <= 0) ? 0 : sobelValues[maxIndex-1];
					y1 = sobelValues[maxIndex];
					y2 = (maxIndex >= stripeLength-3) ? 0 : sobelValues[maxIndex+1];

					//formula for calculating the x-coordinate of the vertex of a parabola, given 3 points with equal distances 
					//(xv means the x value of the vertex, d the distance between the points): 
					//xv = x1 + (d / 2) * (y2 - y0)/(2*y1 - y0 - y2)

					double pos = (y2 - y0) / (4*y1 - 2*y0 - 2*y2 ); //d = 1 because of the normalization and x1 will be added later
					
					// This would be a valid check, too
					//if (std::isinf(pos)) {
					//	// value is infinity
					//	continue;
					//}

					if (pos!=pos) {
						// value is not a number
						continue;
					}

					CvPoint2D32f edgeCenter; //exact point with subpixel accuracy
					int maxIndexShift = maxIndex - (stripeLength>>1);

					//shift the original edgepoint accordingly
					edgeCenter.x = (double)p.x + (((double)maxIndexShift+pos) * stripeVecY.x);
					edgeCenter.y = (double)p.y + (((double)maxIndexShift+pos) * stripeVecY.y);

					CvPoint p_tmp;
					p_tmp.x = (int)edgeCenter.x;
					p_tmp.y = (int)edgeCenter.y;
					cvCircle (iplGrabbed, p_tmp, 1, CV_RGB(0,0,255), -1);

					points[j-1].x = edgeCenter.x;
					points[j-1].y = edgeCenter.y;

					if (isFirstStripe)
					{
						IplImage* iplTmp = cvCreateImage( cvSize(100,300), IPL_DEPTH_8U, 1 );
						cvResize( iplStripe, iplTmp, CV_INTER_NN );
						cvShowImage ( "Exercise 8 - Stripe", iplTmp );//iplStripe );
						cvReleaseImage( &iplTmp );
						isFirstStripe = false;
					}

				} // end of loop over edge points of one edge
				cvReleaseImage ( &iplStripe );

				// we now have the array of exact edge centers stored in "points"
				CvMat mat = cvMat ( 1, 6, CV_32FC2, points);
				cvFitLine ( &mat, CV_DIST_L2, 0, 0.01, 0.01, &lineParams[4*i] );
				// cvFitLine stores the calculated line in lineParams in the following way:
				// vec.x, vec.y, point.x, point.y

				CvPoint p;
				p.x=(int)lineParams[4*i+2] - (int)(50.0*lineParams[4*i+0]);
				p.y=(int)lineParams[4*i+3] - (int)(50.0*lineParams[4*i+1]);

				CvPoint p2;
				p2.x = (int)lineParams[4*i+2] + (int)(50.0*lineParams[4*i+0]);
				p2.y = (int)lineParams[4*i+3] + (int)(50.0*lineParams[4*i+1]);

				cvLine (iplGrabbed, p, p2, CV_RGB(0,255,255), 1, 8, 0);

			} // end of loop over the 4 edges

			// so far we stored the exact line parameters and show the lines in the image
			// now we have to calculate the exact corners
			CvPoint2D32f corners[4];

			for (int i=0; i<4; ++i)
			{
				int j = (i+1)%4;
				double x0,x1,y0,y1,u0,u1,v0,v1;
				x0 = lineParams[4*i+2]; y0 = lineParams[4*i+3];
				x1 = lineParams[4*j+2]; y1 = lineParams[4*j+3];

				u0 = lineParams[4*i+0]; v0 = lineParams[4*i+1];
				u1 = lineParams[4*j+0]; v1 = lineParams[4*j+1];

				// (x|y) = p + s * vec
				// s = Ds / D (see cramer's rule)
				// (x|y) = p + (Ds / D) * vec
				// (x|y) = (p * D / D) + (Ds * vec / D)
				// (x|y) = (p * D + Ds * vec) / D
				// (x|y) = a / c;
				double a =  x1*u0*v1 - y1*u0*u1 - x0*u1*v0 + y0*u0*u1;
				double b = -x0*v0*v1 + y0*u0*v1 + x1*v0*v1 - y1*v0*u1;
				double c =  v1*u0-v0*u1;

				if ( fabs(c) < 0.001 ) //lines parallel?
				{
					std::cout << "lines parallel" << std::endl;
					continue;
				}

				a /= c;
				b /= c;

				//exact corner
				corners[i].x = a; 
				corners[i].y = b;
				CvPoint p;
				p.x = (int)corners[i].x;
				p.y = (int)corners[i].y;

				cvCircle (iplGrabbed, p, 5, CV_RGB(i*60,i*60,0), -1);
			} //finished the calculation of the exact corners

			// resultMatrix_005A made global variable

			CvPoint2D32f targetCorners[4];
			targetCorners[0].x = -0.5; targetCorners[0].y = -0.5;
			targetCorners[1].x =  5.5; targetCorners[1].y = -0.5;
			targetCorners[2].x =  5.5; targetCorners[2].y =  5.5;
			targetCorners[3].x = -0.5; targetCorners[3].y =  5.5;

			//create and calculate the matrix of perspective transform
			CvMat* projMat = cvCreateMat (3, 3, CV_32F );
			cvWarpPerspectiveQMatrix ( corners, targetCorners, projMat);

			//create image for the marker
			CvSize markerSize;
			markerSize.width  = 6;
			markerSize.height = 6;
			IplImage* iplMarker = cvCreateImage( markerSize, IPL_DEPTH_8U, 1 );

			//change the perspective in the marker image using the previously calculated matrix
			cvWarpPerspective(iplConverted, iplMarker, projMat, CV_WARP_FILL_OUTLIERS,  cvScalarAll(0));
			
			cvThreshold(iplMarker, iplMarker, bw_thresh, 255, CV_THRESH_BINARY);

//now we have a B/W image of a supposed Marker

			// check if border is black
			int code = 0;
			for (int i = 0; i < 6; ++i)
			{
				int pixel1 = ((unsigned char*)(iplMarker->imageData + 0*iplMarker->widthStep + i))[0]; //top
				int pixel2 = ((unsigned char*)(iplMarker->imageData + 5*iplMarker->widthStep + i))[0]; //bottom
				int pixel3 = ((unsigned char*)(iplMarker->imageData + i*iplMarker->widthStep))[0]; //left
				int pixel4 = ((unsigned char*)(iplMarker->imageData + i*iplMarker->widthStep + 5))[0]; //right
				if ( ( pixel1 > 0 ) || ( pixel2 > 0 ) || ( pixel3 > 0 ) || ( pixel4 > 0 ) )
				{
					code = -1;
					break;
				}
			}

			if ( code < 0 ) continue;

			//copy the BW values into cP
			int cP[4][4];
			for ( int i=0; i < 4; ++i)
			{
				for ( int j=0; j < 4; ++j)
				{
					cP[i][j] = ((unsigned char*)(iplMarker->imageData + (i+1)*iplMarker->widthStep + (j+1) ))[0];
					cP[i][j] = (cP[i][j]==0) ? 1 : 0; //if black then 1 else 0
				}
			}

			//save the ID of the marker
			int codes[4];
			codes[0] = codes[1] = codes[2] = codes[3] = 0;
			for (int i=0; i < 16; i++)
			{
				int row = i>>2;
				int col = i%4;

				codes[0] <<= 1;
				codes[0] |= cP[row][col]; // 0°

				codes[1] <<= 1;
				codes[1] |= cP[3-col][row]; // 90°

				codes[2] <<= 1;
				codes[2] |= cP[3-row][3-col]; // 180°

				codes[3] <<= 1;
				codes[3] |= cP[col][3-row]; // 270°
			}

			if ( (codes[0]==0) || (codes[0]==0xffff) )
				continue;

			//account for symmetry
			code = codes[0];
			int angle = 0;
			for ( int i=1; i<4; ++i )
			{
				if ( codes[i] < code )
				{
					code = codes[i];
					angle = i;
				}
			}

			//correct order of corners
			if(angle != 0)
			{
				CvPoint2D32f corrected_corners[4];
				for(int i = 0; i < 4; i++)	corrected_corners[(i + angle)%4] = corners[i];
				for(int i = 0; i < 4; i++)	corners[i] = corrected_corners[i];
			}

     		printf ("Found: %04x\n", code);
			if ( isFirstMarker )
			{
				cvShowImage ( "Marker", iplMarker );
				isFirstMarker = false;
			}

			// transfer camera coords to screen coords
			for(int i = 0; i<4; i++)
			{
				corners[i].x -= width/2;
				corners[i].y = -corners[i].y + height/2;
			}

// Added in Exercise 9 - Start *****************************************************************
			//use different display function? define that the marker is smaller?
			if(code == 0x005a)
			{
				estimateSquarePose( resultMatrix_005A, corners, 0.045 );
				detected[0]=true;//display cover
			}
			if(code == 0x0272)
			{
				estimateSquarePose( resultMatrix_0272, corners, 0.045 );
				detected[3]=true;//display image
			}
			if(code == 0x1c44)
			{
				estimateSquarePose( resultMatrix_1c44, corners, 0.045 );
				detected[2]=true;//turn right
			}
			if(code == 0x1228)
			{
				estimateSquarePose( resultMatrix_1228, corners, 0.045 );
				detected[1]=true;//turn left
			}
// Added in Exercise 9 - End *****************************************************************

			/*
				for (int i = 0; i<4; ++i) {
					for (int j = 0; j<4; ++j) {
						cout << setw(6);
						cout << setprecision(4);
						cout << resultMatrix_005A[4*i+j] << " ";
					}
					cout << "\n";
				}
				cout << "\n";
				float x,y,z;
				x = resultMatrix_005A[3];
				y = resultMatrix_005A[7];
				z = resultMatrix_005A[11];
				cout << "length: " << sqrt(x*x+y*y+z*z) << "\n";
				cout << "\n";
			*/

			cvReleaseMat (&projMat);

			delete[] rect;
		} // end of if(result->total == 4)
	} // end of loop over contours

	cvShowImage("Exercise 8 - Original Image", iplGrabbed);
	cvShowImage("Exercise 8 - Converted Image", iplThreshold);

	int key = cvWaitKey (10);
	if (key == 27) exit(0);
// Added in Exercise 9 - Start *****************************************************************
	else if (key == 100) debugmode = !debugmode;
	else if (key == 98) balldebug = !balldebug;
// Added in Exercise 9 - End *****************************************************************

	isFirstStripe = true;

	isFirstMarker = true;

	cvReleaseImage (&iplConverted);
	cvReleaseImage (&iplThreshold);

	cvClearMemStorage ( memStorage );
	glutPostRedisplay();
}

void cleanup() 
{
	cvReleaseMemStorage (&memStorage);

	cvReleaseCapture (&cap);
	cvDestroyWindow ("Exercise 8 - Original Image");
	cvDestroyWindow ("Exercise 8 - Converted Image");
	cvDestroyWindow ("Exercise 8 - Stripe");
	cvDestroyWindow ("Marker");
	cout << "Finished\n";
}

// Added in Exercise 9 - Start *****************************************************************
void drawMarker( int num )
{
	glRotatef( 90, 0, 0, 1 );
	glScalef(0.03, 0.03, 0.03);

	//display texture
	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, texture[num]);	
		glBegin(GL_QUADS);													// Begin Drawing Quads					// (Example Code... Can Be Removed)
		// Front Face																							// (Example Code... Can Be Removed)
		glTexCoord2f(1.0f, 0.0f); glVertex3f(-1.0f, -1.0f,  1.0f);												// (Example Code... Can Be Removed)
		glTexCoord2f(0.0f, 0.0f); glVertex3f( 1.0f, -1.0f,  1.0f);												// (Example Code... Can Be Removed)
		glTexCoord2f(0.0f, 1.0f); glVertex3f( 1.0f,  1.0f,  1.0f);												// (Example Code... Can Be Removed)
		glTexCoord2f(1.0f, 1.0f); glVertex3f(-1.0f,  1.0f,  1.0f);											// (Example Code... Can Be Removed)
	glEnd();
}
// Added in Exercise 9 - End *****************************************************************
void detectFlipping()
{
	bool rightdetected = detected[1];
	bool leftdetected = detected[2]; 
	if (!leftdetected){} //checkleft
	if (!rightdetected) {}//checkright
	if(!leftdetected && rightdetected) turn=1;
	if (!rightdetected && leftdetected) turn=2;
	if (leftdetected)
	{
		if (turn == 1)//flipanimation
		{
			turn = 0;
			if (picture<5)
			{
				picture=picture++;
			}
		}
		for (int x=0; x<4; ++x)
		{
			for (int y=0; y<4; ++y)
			{
				resultTransposedMatrix[x*4+y] = resultMatrix_1c44[y*4+x];
			}
		}
		glLoadMatrixf( resultTransposedMatrix );
		drawMarker(1);
	}
	if (rightdetected)
	{
		if (turn == 2)//flipanimation
		{
			turn = 0;
			if (picture>3)
			{
				picture=picture--;
			}

		}
		for (int x=0; x<4; ++x)
		{
			for (int y=0; y<4; ++y)
			{
				resultTransposedMatrix[x*4+y] = resultMatrix_1228[y*4+x];
			}
		}
		glLoadMatrixf( resultTransposedMatrix );
		drawMarker(2);
	}
}

void display() 
{
    // clear buffers
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    // draw background image
    glDisable( GL_DEPTH_TEST );

    glMatrixMode( GL_PROJECTION );
    glPushMatrix();
    glLoadIdentity();
    gluOrtho2D( 0.0, width, 0.0, height );

    glRasterPos2i( 0, height-1 );
    glDrawPixels( width, height, GL_BGR_EXT, GL_UNSIGNED_BYTE, bkgnd );

    glPopMatrix();

    glEnable(GL_DEPTH_TEST);

    // move to origin
    glMatrixMode( GL_MODELVIEW );
    if (detected[0])
	{
		for (int x=0; x<4; ++x)
		{
			for (int y=0; y<4; ++y)
			{
				resultTransposedMatrix[x*4+y] = resultMatrix_005A[y*4+x];
			}
		}
		glLoadMatrixf( resultTransposedMatrix );
		drawMarker( 0 );
	}
	else if (detected[3])
	{
		for (int x=0; x<4; ++x)
		{
			for (int y=0; y<4; ++y)
			{
				resultTransposedMatrix[x*4+y] = resultMatrix_0272[y*4+x];
			}
		}
		glLoadMatrixf( resultTransposedMatrix );
		drawMarker( picture );
		detectFlipping();
	}
	//draw image


	//rotateToMarker(resultMatrix_005A, resultMatrix_0272, 0x005a);
	//rotateToMarker(resultMatrix_0272, resultMatrix_005A, 0x0272);


	//drawBall
	glLoadIdentity();
	glTranslatef(ballpos.x, ballpos.y + 0.024, ballpos.z);
	//glColor4f(1,0,0,1);
	glutSolidSphere(0.005, 10, 10);
// Added in Exercise 9 - End *****************************************************************

    // redraw
    glutSwapBuffers();
}
//calculates the flipping actions

void resize( int w, int h) 
{
//    width = w;
  //  height = h;

    // set a whole-window viewport
    glViewport( 0, 0, width, height );

    // create a perspective projection matrix
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
	// Note: Just setting the Perspective is an easy hack. In fact, the camera should be calibrated.
	// With such a calibration we would get the projection matrix. This matrix could then be loaded 
	// to GL_PROJECTION.
	// If you are using another camera (which you'll do in most cases), you'll have to adjust the FOV
	// value. How? Fiddle around: Move Marker to edge of display and check if you have to increase or 
	// decrease.
    gluPerspective( camangle, ((double)width/(double)height), 0.01, 100 );

    // invalidate display
    glutPostRedisplay();
}

int main(int argc, char* argv[]) 
{
	cout << "Startup\n";

    // initialize the window system
    glutInit( &argc, argv );
    glutInitWindowSize( width, height );
    glutInitDisplayMode( GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH );
    glutCreateWindow("AR Exercise 8 - Combine");

    // initialize the GL library

    // pixel storage/packing stuff
    glPixelStorei( GL_PACK_ALIGNMENT,   1 );
    glPixelStorei( GL_UNPACK_ALIGNMENT, 1 );
    glPixelZoom( 1.0, -1.0 );

    // enable and set colors
    glEnable( GL_COLOR_MATERIAL );
    glClearColor( 0, 0, 0, 1.0 );

    // enable and set depth parameters
    glEnable( GL_DEPTH_TEST );
    glClearDepth( 1.0 );

    // light parameters
    GLfloat light_pos[] = { 1.0, 1.0, 1.0, 0.0 };
    GLfloat light_amb[] = { 0.2, 0.2, 0.2, 1.0 };
    GLfloat light_dif[] = { 0.7, 0.7, 0.7, 1.0 };

    // enable lighting
    glLightfv( GL_LIGHT0, GL_POSITION, light_pos );
    glLightfv( GL_LIGHT0, GL_AMBIENT,  light_amb );
    glLightfv( GL_LIGHT0, GL_DIFFUSE,  light_dif );
    glEnable( GL_LIGHTING );
    glEnable( GL_LIGHT0 );

    // make functions known to GLUT
    glutDisplayFunc( display );
    glutReshapeFunc( resize  );
    glutIdleFunc( idle );

    // setup OpenCV
    init();

    // for tracker debugging...
    //while (1) idle();

// Added in Exercise 9 - Start *****************************************************************
	std::cout << glGetString(GL_VENDOR) << std::endl;
	std::cout << glGetString(GL_RENDERER) << std::endl;
	std::cout << glGetString(GL_VERSION) << std::endl;
//	std::cout << glGetString(GL_EXTENSIONS) << std::endl;
	std::cout << "Press 'd' for printing debug information. Press 'b' for ball debug information." << std::endl;
// Added in Exercise 9 - End *****************************************************************

    // start the action
    glutMainLoop();
  
    return 0;
}



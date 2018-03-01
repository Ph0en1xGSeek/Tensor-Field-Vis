
#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <string.h>
#include <stdio.h>
#include "f2c.h"
#include <cmath>
#include <gl\glew.h>
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <conio.h>
#define M 59280
#define X 160
#define Y 190
#define Z 148
const double eps = 0.001;

extern"C"
{
#include <clapack.h>
}

using namespace std;

typedef union {
	float f;
	char c[4];
}FLOAT_CONV;

int up_ang = 0;
int left_ang = 0;
float zoom = 0.6f;
bool isMouseDown = false;
int clickx, clicky;
///---------------------------------------------------------------
float arr[9];
doublereal A[9];
doublereal* wr[X][Y][Z];
doublereal* vl[X][Y][Z];
double cl[X][Y][Z];
double cp[X][Y][Z];
double cs[X][Y][Z];
double fa[X][Y][Z];
bool flags[X][Y][Z];
double mxcl = 0, mxcp = 0, mxcs = 0, mxfa = 0;
double micl = 1, micp = 1, mics = 1, mifa = 1;
///--------------------------------------------------------------------
float mtrl_ambient[] = { 0.19225, 0.19225, 0.19225, 1.0 };  // 僔儖僶乕(攚宨1)
float mtrl_diffuse[] = { 0.50754, 0.50754, 0.50754, 1.0 };
float mtrl_specular[] = { 0.21, 0.21, 0.21, 1.0 };
float mtrl_shininess[] = { 0.4 };
float no_mat[] = { 1.0,1.0,1.0,0 };
float zero_mat[] = { 0,0,0,1 };
float one_mat[] = { 1.0,1.0,1.0,1 };
int color_mode = 1;
///====================================================================
int transX = 0, transY = 0;


float __ltobf(float data)
{
	FLOAT_CONV d1, d2;

	d1.f = data;

	d2.c[0] = d1.c[3];
	d2.c[1] = d1.c[2];
	d2.c[2] = d1.c[1];
	d2.c[3] = d1.c[0];
	return d2.f;
}

void init(void)
{
	glClearColor(0, 0, 0, 0);
	glShadeModel(GL_SMOOTH);

	float light0_position[] = { 500, 500,  500,  -1 };
	float light0_ambient[] = { 0.7, 0.7, 0.7, 1.0 };
	float light0_diffuse[] = { 0.7, 0.7, 0.7, 1.0 };
	float light0_specular[] = { 0.3, 0.3, 0.3, 1.0 };

	float light1_position[] = { -500, 500,  500,  -1 };
	float light1_ambient[] = { 0.3, 0.3, 0.3, 1.0 };
	float light1_diffuse[] = { 0.3, 0.3, 0.3, 1.0 };
	float light1_specular[] = { 1.0, 1.0, 1.0, 0.0 };


	glLightfv(GL_LIGHT0, GL_POSITION, light0_position);
	glLightfv(GL_LIGHT0, GL_AMBIENT, light0_ambient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light0_diffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light0_specular);

	glLightfv(GL_LIGHT1, GL_POSITION, light1_position);
	glLightfv(GL_LIGHT1, GL_AMBIENT, light1_ambient);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, light1_diffuse);
	glLightfv(GL_LIGHT1, GL_SPECULAR, light1_specular);


	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);
	//glDisable(GL_COLOR_MATERIAL);
	glEnable(GL_BLEND);
	glMateriali(GL_FRONT, GL_SHININESS, 100);
	glEnable(GL_NORMALIZE);
	//glEnable(GL_POINT_SMOOTH);
	//glEnable(GL_LINE_SMOOTH);
	//glEnable(GL_POLYGON_SMOOTH);
	//glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
	//glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
	//glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
	glEnable(GL_MULTISAMPLE_ARB);//多重采样
}

int j = 70;

void display(void)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(0.0, 0.0, 0.0, 1.0);
	glLoadIdentity();
	//gluLookAt(0, Y, 0, 0, 0.0, 0, -1, 0, 0);


	glPushMatrix();
	glRotatef(up_ang + 90, 1, 0, 0);
	glRotatef(left_ang, 0, 1, 0);
	glScalef(zoom, zoom, zoom);
	glTranslatef(-(Z*1.9*1.6) / 2, 0, 0);
	//glTranslatef(0, -Y*1.6*1.48 / 2, 0);
	glTranslatef(0, 0, -X*1.48*1.9 / 2);
	glTranslated(transX, 0, transY);
	for (int i = 0; i < X; i++)
	{
		//for (int j = 0; j < Y; j++)
		//{
		for (int k = 0; k < Z; k++)
		{
			if (flags[i][j][k] == false) continue;
			double a = wr[i][j][k][0], b = wr[i][j][k][1], c = wr[i][j][k][2];

			double mx = -999.999;
			//standarize
			mx = max(a, b);
			mx = max(mx, c);
			a /= mx;
			b /= mx;
			c /= mx;
			no_mat[0] = cl[i][j][k];
			no_mat[1] = cl[i][j][k];
			no_mat[2] = cl[i][j][k];
			no_mat[3] = 1;
			glMaterialfv(GL_FRONT, GL_AMBIENT, no_mat);
			glMaterialfv(GL_FRONT, GL_DIFFUSE, no_mat);
			glMaterialfv(GL_FRONT, GL_SPECULAR, no_mat);
			glMaterialfv(GL_FRONT, GL_SHININESS, zero_mat);
			//glColor3f(cl[i][j][k],cl[i][j][k],cl[i][j][k]);
			glBegin(GL_QUADS);
			glVertex3f(k*1.9*1.6 + 1.9*1.6 / 2, 1.5, i*1.9*1.48 + 1.9*1.48 / 2);
			glVertex3f(k*1.9*1.6 - 1.9*1.6 / 2, 1.5, i*1.9*1.48 + 1.9*1.48 / 2);
			glVertex3f(k*1.9*1.6 - 1.9*1.6 / 2, 1.5, i*1.9*1.48 - 1.9*1.48 / 2);
			glVertex3f(k*1.9*1.6 + 1.9*1.6 / 2, 1.5, i*1.9*1.48 - 1.9*1.48 / 2);
			glEnd();
			if (fa[i][j][k] < 0.42)
			{
				continue;
			}
			double tmp[9];
			double rotate[16];
			for (int ii = 0; ii < 9; ii++)
			{
				tmp[ii] = vl[i][j][k][ii];
				//把每个特征向量模变为1
				if (ii % 3 == 2)
				{
					double tot = sqrt(tmp[ii - 2] * tmp[ii - 2] + tmp[ii - 1] * tmp[ii - 1] + tmp[ii] * tmp[ii]);
					tmp[ii - 2] /= (tot);
					tmp[ii - 1] /= (tot);
					tmp[ii] /= (tot);
				}
			}
			//========================
			rotate[0] = tmp[0];
			rotate[1] = tmp[1];
			rotate[2] = tmp[2];
			rotate[4] = tmp[3];
			rotate[5] = tmp[4];
			rotate[6] = tmp[5];
			rotate[8] = tmp[6];
			rotate[9] = tmp[7];
			rotate[10] = tmp[8];
			rotate[3] = rotate[7] = rotate[11] = rotate[12] = rotate[13] = rotate[14] = 0.0;
			rotate[15] = 1.0;
			//=======================
			double red = (cl[i][j][k] - micl) / (mxcl - micl);
			double green = (cp[i][j][k] - micp) / (mxcp - micp);
			double blue = (cs[i][j][k] - mics) / (mxcs - mics);
			double alpha = (fa[i][j][k] - mifa) / (mxfa - mifa);
			//---------------------------------------------------------------------配色方案1------------------------------
			if (color_mode == 1)
			{
				mtrl_ambient[0] = cl[i][j][k], mtrl_ambient[1] = cp[i][j][k], mtrl_ambient[2] = cs[i][j][k], mtrl_ambient[3] = 1.0;
				mtrl_diffuse[0] = cl[i][j][k], mtrl_diffuse[1] = cp[i][j][k], mtrl_diffuse[2] = cs[i][j][k], mtrl_diffuse[3] = 1.0;
				mtrl_specular[0] = cl[i][j][k], mtrl_specular[1] = cp[i][j][k], mtrl_specular[2] = cs[i][j][k], mtrl_specular[3] = 1.0;
			}
			//----------------------------------------------------------------------------------配色方案2-----------------------------------------
			else if (color_mode == 2)
			{
				mtrl_ambient[0] = red, mtrl_ambient[1] = green, mtrl_ambient[2] = blue, mtrl_ambient[3] = 1.0;
				mtrl_diffuse[0] = red, mtrl_diffuse[1] = green, mtrl_diffuse[2] = blue, mtrl_diffuse[3] = 1.0;
				mtrl_specular[0] = red, mtrl_specular[1] = green, mtrl_specular[2] = blue, mtrl_specular[3] = 1.0;

			}
			//-------------------------------------------------------------------------------------配色方案3--------------------------------------
			else if (color_mode == 3)
			{

				if (cl[i][j][k] >= cp[i][j][k] && cl[i][j][k] >= cs[i][j][k])
				{
					//glColor4d(red, 0.1, 0.1, 1);
					mtrl_ambient[0] = red, mtrl_ambient[1] = 0.1, mtrl_ambient[2] = 0.1, mtrl_ambient[3] = 1.0;
					mtrl_diffuse[0] = red, mtrl_diffuse[1] = 0.1, mtrl_diffuse[2] = 0.1, mtrl_diffuse[3] = 1.0;
					mtrl_specular[0] = red, mtrl_specular[1] = 0.1, mtrl_specular[2] = 0.1, mtrl_specular[3] = 1.0;
				}
				else if (cp[i][j][k] >= cl[i][j][k] && cp[i][j][k] >= cs[i][j][k])
				{
					//glColor4d(0.1, green, 0.1, 1);
					mtrl_ambient[0] = 0.1, mtrl_ambient[1] = green, mtrl_ambient[2] = 0.1, mtrl_ambient[3] = 1.0;
					mtrl_diffuse[0] = 0.1, mtrl_diffuse[1] = green, mtrl_diffuse[2] = 0.1, mtrl_diffuse[3] = 1.0;
					mtrl_specular[0] = 0.1, mtrl_specular[1] = green, mtrl_specular[2] = 0.1, mtrl_specular[3] = 1.0;
				}
				else
				{
					//glColor4d(0.1, 0.1, blue, 1);
					mtrl_ambient[0] = 0.1, mtrl_ambient[1] = 0.1, mtrl_ambient[2] = blue, mtrl_ambient[3] = 1.0;
					mtrl_diffuse[0] = 0.1, mtrl_diffuse[1] = 0.1, mtrl_diffuse[2] = blue, mtrl_diffuse[3] = 1.0;
					mtrl_specular[0] = 0.1, mtrl_specular[1] = 0.1, mtrl_specular[2] = blue, mtrl_specular[3] = 1.0;
				}
			}
			//----------------------------------------------------------------------------------------------------------------------------------
			glMaterialfv(GL_FRONT, GL_AMBIENT, mtrl_ambient);
			glMaterialfv(GL_FRONT, GL_DIFFUSE, mtrl_diffuse);
			glMaterialfv(GL_FRONT, GL_SPECULAR, one_mat);
			glMaterialfv(GL_FRONT, GL_SHININESS, one_mat);



			//---------------------------------新方法
			glPushMatrix();
			glTranslated(k*1.9*1.6, -1.5, i*1.9*1.48);
			glMultMatrixd(rotate);
			glScaled(a,b,c);
			glutSolidSphere(3.0, 30, 30);
			//glutSolidCube(4);
			glPopMatrix();
			//--------------------------------
		}
		//}
	}
	glPopMatrix();
	glFlush();
	glutSwapBuffers();
	//glutPostRedisplay();
}

void reshape(int w, int h)
{
	glViewport(0, 0, (GLsizei)w, (GLsizei)h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	//gluPerspective(90, w / h, 0, 60);
	glOrtho(-200.0, 200.0, -200.0, 200.0, 200.0, -200.0);
	glMatrixMode(GL_MODELVIEW);
	//glLoadIdentity();
}

void keyboard(unsigned char key, int x, int y)
{
	if (key == 'w')
	{
		up_ang++;
	}
	else if (key == 's')
	{
		up_ang--;
	}
	else if (key == 'a')
	{
		left_ang++;
	}
	else if (key == 'd')
	{
		left_ang--;
	}
	else if (key == '+')
	{
		printf("\b\b\b");
		j++;
		j = min(j, Y - 1);
		printf("%3d", j);
	}
	else if (key == '-')
	{
		printf("\b\b\b");
		j--;
		j = max(j, 0);
		printf("%3d", j);
	}
	else if (key == '1')
	{
		color_mode = 1;
	}
	else if (key == '2')
	{
		color_mode = 2;
	}
	else if (key == '3')
	{
		color_mode = 3;
	}
	//glutPostRedisplay();
}

void mouse(int key, int state, int x, int y)
{
	if (state == GLUT_UP && key == GLUT_WHEEL_UP)
	{
		zoom += 0.03;
	}
	if (state == GLUT_UP && key == GLUT_WHEEL_DOWN)
	{
		zoom -= 0.03;
	}
	zoom = max(zoom, 0);
	//glutPostRedisplay();
	if (key == GLUT_LEFT_BUTTON && state == GLUT_DOWN)
	{
		isMouseDown = true;
		clickx = x;
		clicky = y;
	}
	else if (key == GLUT_LEFT_BUTTON && state == GLUT_UP)
	{
		isMouseDown = false;
	}
}

void motion(int x, int y)
{
	if (isMouseDown)
	{
		if (y > clicky)
		{
			transY -= (clicky - y);
			clicky = y;
			//printf("r_ang %lf\n", r_angle);
			//            printf("up\n");
		}
		else if (y < clicky)
		{
			transY -= (clicky - y);
			clicky = y;
			//printf("r_ang %lf\n", r_angle);
			//            printf("down\n");
		}

		if (x < clickx)
		{
			transX -= (clickx - x);
			clickx = x;
			//printf("s_ang %lf\n", s_angle);
			//            printf("up\n");
		}
		else if (x > clickx)
		{
			transX -= (clickx - x);
			clickx = x;
			//printf("s_ang %lf\n", s_angle);
			//            printf("down\n");
		}
	}
}

void myIdle(void)
{
	glutPostRedisplay();
}


void specialKeys(int key, int x, int y)
{
	switch (key)
	{
	case GLUT_KEY_UP:
		transY--;
		break;
	case GLUT_KEY_DOWN:
		transY++;
		break;
	case GLUT_KEY_LEFT:
		transX--;
		break;
	case GLUT_KEY_RIGHT:
		transX++;
		break;
	}
}

int main(int argc, char** argv)
{
	///------------------------------------------------
	memset(flags, false, sizeof(flags));
	integer info;
	float flag;
	char jobvl = 'V';
	char jobvr = 'V';
	integer n = 3;
	doublereal *a = (doublereal*)A;
	integer lda = 3;

	doublereal wi[3];

	integer ldvr = 3;
	doublereal* vr = new doublereal[n * ldvr];

	integer ldvl = 3;

	integer lwork = n * 4;
	doublereal *work = new doublereal[lwork];
	int cnt = 0;
	ifstream fin("gk2-rcc-mask.raw", ios::binary);
	if (!fin)
	{
		cout << "error" << endl;
	}
	else
	{
		for (int i = 0; i < X; i++)
		{
			for (int j = 0; j < Y; j++)
			{
				for (int k = 0; k < Z; k++)
				{
					fin.read((char*)&flag, sizeof(float));
					fin.read((char*)&arr[0], sizeof(float));
					fin.read((char*)&arr[1], sizeof(float));
					fin.read((char*)&arr[2], sizeof(float));
					fin.read((char*)&arr[4], sizeof(float));
					fin.read((char*)&arr[5], sizeof(float));
					fin.read((char*)&arr[8], sizeof(float));
					flag = __ltobf(flag);
					if (fabs(flag - 0.0f) < eps)
					{
						//                        cout << "ignore " << cnt++ << endl;
						flags[i][j][k] = false;
					}
					else
					{
						flags[i][j][k] = true;
						vl[i][j][k] = new doublereal[n * ldvl];
						wr[i][j][k] = new doublereal[3];
						A[0] = (double)__ltobf(arr[0]);
						A[1] = (double)__ltobf(arr[1]);
						A[2] = (double)__ltobf(arr[2]);
						A[4] = (double)__ltobf(arr[4]);
						A[5] = (double)__ltobf(arr[5]);
						A[8] = (double)__ltobf(arr[8]);
						A[3] = (double)__ltobf(arr[1]);
						A[6] = (double)__ltobf(arr[2]);
						A[7] = (double)__ltobf(arr[5]);
						for (int ii = 0; ii < 9; ii++)//-0.00000 != 0.00000
						{
							if (fabs(A[ii]) < 0.00001)
								A[ii] = 0.0000;
						}
						//                            puts("\nmatrix:");
						//                        for (int ii = 0; ii < 3; ii++)
						//                        {
						//                        for (int jj = 0; jj < 3; jj++)
						//                        {
						//                        printf("%10.5lf", A[ii * n + jj]);
						//                        }
						//                        printf("\n");
						//                        }
						dgeev_(&jobvl, &jobvr, &n, a, &lda, wr[i][j][k], wi, vl[i][j][k], &ldvl, vr, &ldvr, work, &lwork, &info);
						//                        printf("%d: \n", cnt++);
						//                        printf("Info: %d\n", info);
						//                        printf("eigenvalue:\n");
						//                        for (int ii = 0; ii < 3; ii++)
						//                        printf("%lf  ", wr[i][j][k][ii]);
						//                        printf("\n");

						//                        for (int ii = 0; ii < 3; ii++)
						//                        {
						//                        printf("eigenvector: ");
						//                        for (int jj = 0; jj < 3; jj++)
						//                        {
						//                        printf("%lf ", vl[i][j][k][ii * n + jj]);
						//                        }
						//                        printf("\n");
						//                        }
						/*for (int ii = 0; ii < 2; ii++)
						{
						doublereal ans = .000000;
						for (int jj = 0; jj < 3; jj++)
						{
						ans += vl[i][j][k][ii * n + jj] * vl[i][j][k][(ii+ 1) * n + jj];
						}
						if (fabs(ans) > eps)
						printf("%d %d %d Error: %lf\n",i, j, k, ans);
						}
						//Check the data*/
						/*printf("unreal : %d %d %e \n", wi[0], wi[1], wi[2]);*/



						double a = wr[i][j][k][0], b = wr[i][j][k][1], c = wr[i][j][k][2];
						if (a < 0 || b < 0 || c < 0)
						{
							flags[i][j][k] = false;
							continue;
						}
						if (a < b)
						{
							swap(a, b);
						}
						if (a < c)
						{
							swap(a, c);
						}
						if (b < c)
						{
							swap(b, c);
						}
						//                        if(a < 0 || b < 0 || c < 0)
						//                        cout << a << ' ' << b << ' ' << c << endl << endl;
						///Color
						cl[i][j][k] = (a - b) / (a + b + c);
						cp[i][j][k] = 2.0 * (b - c) / (a + b + c);
						cs[i][j][k] = 3.0 * c / (a + b + c);
						double avg = (a + b + c) / 3.0;
						fa[i][j][k] = sqrt(3 * ((a - avg)*(a - avg) + (b - avg)*(b - avg) + (c - avg)*(c - avg))) / sqrt(2 * (a*a + b*b + c*c));
						if (fa[i][j][k] >= 0.42)
						{
							mxcl = max(mxcl, cl[i][j][k]);
							micl = min(micl, cl[i][j][k]);
							mxcp = max(mxcp, cp[i][j][k]);
							micp = min(micp, cp[i][j][k]);
							mxcs = max(mxcs, cs[i][j][k]);
							mics = min(mics, cs[i][j][k]);
							mxfa = max(mxfa, fa[i][j][k]);
							mifa = min(mifa, fa[i][j][k]);
						}
					}

				}
			}
		}
	}
	printf("w/a/s/d：旋转切片\n");
	printf("↑/↓/←/→：移动切片\n");
	printf("鼠标拖动：移动切片\n");
	printf("鼠标滚轮：缩放切片\n");
	printf("+/-：移动截面\n\n\n");
	printf("当前切片：%3d", j);
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH | GLUT_MULTISAMPLE);
	glutInitWindowSize(500, 500);
	glutInitWindowPosition(100, 100);
	glutCreateWindow("Helix");
	init();
	glutDisplayFunc(display);
	glutReshapeFunc(reshape);
	glutKeyboardFunc(keyboard);
	glutSpecialFunc(specialKeys);
	glutMouseFunc(mouse);
	glutMotionFunc(motion);
	glutIdleFunc(myIdle);
	glutMainLoop();
	return 0;
}

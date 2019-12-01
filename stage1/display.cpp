#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <ctime>

#include "constants.hpp"
#include "simulate.hpp"

simulation *sim;
 
void reshape(int width, int height){
    glViewport(0,0,width,height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-WIDTH/2,WIDTH/2-1,-HEIGHT/2,HEIGHT/2-1,-1,1);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}
 
void init(void){
    glClearColor(0.0,0.0,0.0,1.0);
    glPointSize(2.0);
}
 
void Timer(int ex)
{
    glutPostRedisplay();
    glutTimerFunc(1,Timer,0);
}

void DrawCircle(float cx, float cy, float r, int num_segments)
{
    glBegin(GL_LINE_LOOP);
    for (int ii = 0; ii < num_segments; ii++)   {
        float theta = 2.0f * 3.1415926f * float(ii) / float(num_segments);//get the current angle 
        float x = r * cosf(theta);//calculate the x component 
        float y = r * sinf(theta);//calculate the y component 
        glVertex2f(x + cx, y + cy);//output vertex 
    }
    glEnd();
}

void DrawWals()
{
    glBegin(GL_LINE_LOOP);
    glVertex2f(LEFT_WALL, BOTTOM_WALL);
    glVertex2f(LEFT_WALL, TOP_WALL);
    glEnd();

    glBegin(GL_LINE_LOOP);
    glVertex2f(RIGHT_WALL, BOTTOM_WALL);
    glVertex2f(RIGHT_WALL, TOP_WALL);
    glEnd();

    glBegin(GL_LINE_LOOP);
    glVertex2f(LEFT_WALL, BOTTOM_WALL);
    glVertex2f(RIGHT_WALL, BOTTOM_WALL);
    glEnd();

    glBegin(GL_LINE_LOOP);
    glVertex2f(LEFT_WALL, TOP_WALL);
    glVertex2f(RIGHT_WALL, TOP_WALL);
    glEnd();
}
  
void display(void)
{
    static int counter = 0;

    sim->step();

    counter = (counter + 1) % DRAW_STEP;

    if (counter != 0)
        return;

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glColor4f(0.0,1.0,1.0,1.0);

    DrawWals();

    for (int i = 0; i < sim->balls().size(); i++)
        DrawCircle(sim->balls()[i].position.x, sim->balls()[i].position.y, sim->balls()[i].radius(), 50);

    glutSwapBuffers();
}
 
void idle(void){
/* do nothing */
}
 
int main(int argc, char **argv){
    simulation s(1000, time(NULL));
    sim = &s;

    glutInit(&argc,argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
    glutInitWindowPosition(0,0);
    glutInitWindowSize(WIDTH,HEIGHT);
    glutCreateWindow(argv[0]);
    init();
    glutIdleFunc(idle);
    glutReshapeFunc(reshape);
    glutDisplayFunc(display);
    glutTimerFunc(0,Timer,0);
    glutMainLoop();
    return(1);
}

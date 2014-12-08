#include<iostream>
#include<math.h>
#include "hw4.h"
#include<vector>

stack STACK;
z_buffer ZBUFFER;
matrix * normalized_matrix = new matrix();


void init_stack(){
	STACK.mptr = new matrix[10];
	STACK.head = 0;
}

void init_zbuffer(int x, int y){     //x rownum, y colnum
	ZBUFFER.ele = new double*[x];
	for (int i = 0 ; i < x ; i++){
		ZBUFFER.ele[i] = new double[y];
	}
	for(int i = 0; i < x ; i++){
		for(int j = 0; j< y ; j++){
			ZBUFFER.ele[i][j] = -10; //initializ value to large negative number
		}
	}
}

void copy_matrix(matrix *ptr1, matrix *ptr2){
	for(int i = 0; i < 4 ; i++){
		for(int j = 0; j< 4 ; j++){
			ptr2->ele[i][j] = ptr1->ele[i][j];
		}
	}

}
void push_stack(){
	STACK.head++;
	copy_matrix(&STACK.mptr[STACK.head-1],&STACK.mptr[STACK.head]);		
}

void pop_statck(){
	STACK.head--;
}

matrix * peek_stack(){
	return &STACK.mptr[STACK.head];
}

std::vector<double> CalculateIntersection(linequation lines[4], int y){
    std::vector<double> intersections;
    for(int i = 0; i < 4; i++){
        double temp = (-lines[i].c-lines[i].b * y) / lines[i].a ;
        if(temp <= lines[i].x2 && temp >= lines[i].x1){
            intersections.push_back(temp);
        }
    }
    return intersections;
}

void matrix_multiply(matrix *ptr1, matrix *ptr2, matrix *result){
	for(int i = 0; i < 4 ; i++){
		for(int j = 0; j< 4 ; j++){
			double sum = 0;
			for (int k = 0; k < 4 ; k++){
				sum += ptr1->ele[i][k] * ptr2->ele[k][j];
			}
			result->ele[i][j] = sum;
		}
	}	
}

void vector_matrix_mutiply(point3d *ptr1, matrix * ptr2, point3d *result){
	for(int i = 0 ; i < 4 ; i++){
		double sum = 0;
		for (int j = 0 ; j < 4 ; j++){
			sum+= ptr2->ele[i][j] * ptr1->ele[j];
		}
		result->ele[i] = sum;
	}
}

void multiply_top_stack(matrix * m){
	matrix *copy = new matrix();
	matrix *ptr = peek_stack();
	copy_matrix(ptr,copy);	
	matrix_multiply(copy,m,ptr);
	delete(copy);
}

double get_zbuffer(point2d *ptr){
	int x = ptr->ele[0];
	int y = ptr->ele[1];
	return ZBUFFER.ele[x][y];	
}

void update_zbuffer(point2d *ptr, double z){
	int x = ptr->ele[0];
	int y = ptr->ele[1];
	ZBUFFER.ele[x][y] = z;
}

int compare_update_zbuffer(point2d *ptr, double new_z){
	double old_z = get_zbuffer(ptr);
	if (new_z > old_z){
		update_zbuffer(ptr, new_z);
		return 1;
	}
	else{
		return 0;
	}
}

double deg_to_rad(double deg){
    double pie = 3.1415926;
    return deg*pie/180;
}

matrix * translate(double x, double y, double z){
    matrix *T = new matrix();
    for(int i=0; i<4; i++){
        for(int j=0; j<4; j++){
            T->ele[i][j] = 0;
        }
    }
    T->ele[0][0] = 1;
    T->ele[1][1] = 1;
    T->ele[2][2] = 1;
    T->ele[3][3] = 1;
    T->ele[0][3] = x;
    T->ele[1][3] = y;
    T->ele[2][3] = z;
    return T;
}

void gl_translate(double x, double y, double z){
    matrix * T;
    T = translate(x,y,z);
	//multipy
	multiply_top_stack(T);
	delete(T);
}

void rotate_z(double deg){ // rotate z direction
	matrix *R = new matrix();
    double theta = deg_to_rad(deg);
	//init T according to input
    for(int i=0; i<4; i++){
        for(int j=0; j<4; j++){
            R->ele[i][j] = 0;
        }
    }
    R->ele[0][0] = cos(theta);
    R->ele[1][1] = cos(theta);
    R->ele[0][1] = -sin(theta);
    R->ele[1][0] = sin(theta);
    R->ele[2][2] = 1;
    R->ele[3][3] = 1;
	//multipy 
	multiply_top_stack(R);
	delete(R);
}

void rotate_x(double theta){
	matrix *R = new matrix();
    //init T according to input
    for(int i=0; i<4; i++){
        for(int j=0; j<4; j++){
            R->ele[i][j] = 0;
        }
    }
    R->ele[2][2] = cos(theta);
    R->ele[1][1] = cos(theta);
    R->ele[2][1] = sin(theta);
    R->ele[1][2] = -sin(theta);
    R->ele[0][0] = 1;
    R->ele[3][3] = 1;
	//multipy 
	multiply_top_stack(R);
	delete(R);
}

void rotate_y(double theta){
    matrix *R = new matrix();
    //init T according to input
    for(int i=0; i<4; i++){
        for(int j=0; j<4; j++){
            R->ele[i][j] = 0;
        }
    }
    R->ele[0][0] = cos(theta);
    R->ele[2][2] = cos(theta);
    R->ele[0][2] = sin(theta);
    R->ele[2][0] = -sin(theta);
    R->ele[1][1] = 1;
    R->ele[3][3] = 1;
    //multipy
    multiply_top_stack(R);
    delete(R);
}

matrix * scale(double x, double y, double z){
    matrix *S = new matrix();
    //init T according to input
    for(int i=0; i<4; i++){
        for(int j=0; j<4; j++){
            S->ele[i][j] = 0;
        }
    }
    S->ele[0][0] = x;
    S->ele[1][1] = y;
    S->ele[2][2] = z;
    S->ele[3][3] = 1;
    return S;
}

void gl_scale(double x, double y, double z){
    matrix * S;
    S = scale(x, y, z);
	//multipy
	multiply_top_stack(S);
	delete(S);
}

linequation * returnLineEquation(point2d *ptr1, point2d *ptr2){   // can assume all integer???
    linequation * equation = new linequation();
    equation->a = ptr2->ele[1] - ptr1->ele[1];
    equation->b = ptr1->ele[0] - ptr2->ele[0];
    equation->c = ptr1->ele[1] * (ptr2->ele[0] - ptr1->ele[0]) - ptr1->ele[0] * (ptr2->ele[1] - ptr1->ele[1]);
    if(ptr1->ele[0] < ptr2->ele[0]){//x1<x2
        equation->x1 = ptr1->ele[0];
        equation->x2 = ptr2->ele[0];
    }
    else{
        equation->x1 = ptr2->ele[0];
        equation->x2 = ptr1->ele[0];
    }
    if(ptr1->ele[1] < ptr2->ele[1]){//y1<y2
        equation->y1 = ptr1->ele[1];
        equation->y2 = ptr2->ele[1];
    }
    else{
        equation->y1 = ptr2->ele[1];
        equation->y2 = ptr1->ele[1];
    }
    return equation;
}

planequation * returnPlaneEquation(point3d *ptr1, point3d *ptr2, point3d *ptr3){
    planequation * equation = new planequation();
    equation->a = (ptr2->ele[1] - ptr1->ele[1]) * (ptr3->ele[2] - ptr1->ele[2]) - (ptr3->ele[1] - ptr1->ele[1]) * (ptr2->ele[2] - ptr1->ele[2]);
    equation->b = (ptr2->ele[2] - ptr1->ele[2]) * (ptr3->ele[0] - ptr1->ele[0]) - (ptr3->ele[2] - ptr1->ele[2]) * (ptr2->ele[0] - ptr1->ele[0]);
    equation->c = (ptr2->ele[0] - ptr1->ele[0]) * (ptr3->ele[1] - ptr1->ele[1]) - (ptr3->ele[0] - ptr1->ele[0]) * (ptr2->ele[1] - ptr1->ele[1]);
    equation->d = -(equation->a * ptr1->ele[0] + equation->b * ptr1->ele[1] + equation->c * ptr1->ele[2]);
    return equation;
}

point2d * projection(point3d *ptr){
    point2d * point = new point2d();
    point->ele[0] = int (ptr->ele[0] / ptr->ele[2] + 0.5);  //round double to integer
    point->ele[1] = int (ptr->ele[1] / ptr->ele[2] + 0.5);  //round double to integer
    return point;
}

double getDepth(point2d* ptr, planequation * planeptr){
    double depth;
    depth = (1 / (planeptr->a * ptr->ele[0] + planeptr->b * ptr->ele[1] + planeptr->c))- planeptr->d;
    return depth;
}

void crossProdcut (point3d * x, point3d * y, point3d * result){
    result->ele[3] = 1;
    result->ele[0] = x->ele[1] * y->ele[2] - y->ele[1] * x->ele[2];
    result->ele[1] = y->ele[0] * x->ele[2] - x->ele[0] * y->ele[2];
    result->ele[2] = x->ele[0] * y->ele[1] - y->ele[0] * x->ele[1];
}

void normalization(point3d * eye, point3d * center, point3d * up, double near, double far, double left, double right, double bottom, double top){
    matrix * t;
    matrix * r = new matrix();
    matrix * shear = new matrix();
    matrix * scaleXY = new matrix();
    matrix * scaleZ = new matrix();

    t = translate(-eye->ele[0], -eye->ele[1], -eye->ele[2]);
    double Wmagnitude = sqrt (pow((eye->ele[0] - center->ele[0]),2) + pow((eye->ele[1] - center->ele[1]),2) + pow((eye->ele[2] - center->ele[2]),2) );

    point3d * w = new point3d;
    w->ele[3] = 1;
    w->ele[2] = (eye->ele[2] - center->ele[2]) / Wmagnitude;
    w->ele[1] = (eye->ele[1] - center->ele[1]) / Wmagnitude;
    w->ele[0] = (eye->ele[0] - center->ele[0]) / Wmagnitude;
    
    point3d * u = new point3d;
    crossProdcut(up, w, u);
    double Umagnitude = sqrt (pow(u->ele[0],2) + pow(u->ele[1],2) + pow(u->ele[2],2) );
    u->ele[3] = 1;
    u->ele[2] = u->ele[2] / Umagnitude;
    u->ele[1] = u->ele[1] / Umagnitude;
    u->ele[0] = u->ele[0] / Umagnitude;
    
    point3d * v = new point3d;
    crossProdcut(w, u, v);
    double Vmagnitude = sqrt (pow(v->ele[0],2) + pow(v->ele[1],2) + pow(v->ele[2],2) );
    v->ele[3] = 1;
    v->ele[2] = v->ele[2] / Vmagnitude;
    v->ele[1] = v->ele[1] / Vmagnitude;
    v->ele[0] = v->ele[0] / Vmagnitude;
    
    for(int i=0; i<4; i++){
        for(int j=0; j<4; j++){
            r->ele[i][j] = 0;
        }
    }
    r->ele[3][3] = 1;
    r->ele[0][0] = u->ele[0];
    r->ele[0][1] = u->ele[1];
    r->ele[0][2] = u->ele[2];
    r->ele[1][0] = v->ele[0];
    r->ele[1][1] = v->ele[1];
    r->ele[1][2] = v->ele[2];
    r->ele[2][0] = w->ele[0];
    r->ele[2][1] = w->ele[1];
    r->ele[2][2] = w->ele[2];

    for(int i=0; i<4; i++){
        for(int j=0; j<4; j++){
            shear->ele[i][j] = 0;
        }
    }
    shear->ele[0][0] = 1;
    shear->ele[1][1] = 1;
    shear->ele[2][2] = 1;
    shear->ele[3][3] = 1;
    shear->ele[0][2] = (left + right) / (2 * near);
    shear->ele[1][2] = (top + bottom) / (2 * near);

    scaleXY = scale(((2 * near) / (right - left)), (2 * near) / (top - bottom), 1);
    scaleZ = scale(1 / far, 1 / far, 1 / far);
    
    matrix * temp = new matrix();
    matrix_multiply (t, r, temp);
    matrix_multiply (temp, shear,normalized_matrix);
    matrix_multiply (normalized_matrix, scaleXY, temp);
    matrix_multiply (temp, scaleZ, normalized_matrix);
    delete(temp);
    delete(t);
    delete(r);
    delete(shear);
    delete(scaleXY);
    delete(scaleZ); // need check!!!
    
    matrix * topstack = peek_stack();
    copy_matrix(normalized_matrix, topstack);  // push into the stack??  delete normalized_matrix?
}

void scan(linequation lines[4]){
    int ymin = lines[0].y1, ymax = lines[0].y2;
    
    for(int i = 1; i < 4; i++){
        if(ymin > lines[i].y1)
            ymin = lines[i].y1;
    }
    for(int i = 1; i < 4; i++){
        if(ymax<lines[i].y2){
            ymax = lines[i].y2;
        }
    }
    for(int scan_y = ymin; scan_y <= ymax; scan_y ++){
        std::<vetcor> double intersections = CalculateIntersection(lines, scan_y);
        drawLine(int (intersections.at(0)+0.5), int (intersections.at(1)+0.5), scan_y, scan_y, 0); //round to integer???
    }
}

void drawLine(int x1, int x2, int y1, int y2, int color){   //color: 0 -> black, 1 -> white
    
}

int main(){
	return 0;
}

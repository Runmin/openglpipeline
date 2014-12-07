#include<iostream>
#include "hw4.h"

stack STACK;
z_buffer ZBUFFER;

void init_stack(){
	STACK.mptr = new matrix[10];
	STACK.head = 0;
}

void init_zbuffer(int x, int y){
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

void translate(double x, double y, double z){
	matrix *T = new matrix();
	//init T according to input
	//multipy 
	multiply_top_stack(T);
	delete(T);
}

void rotate_z(double theta){ // rotate z direction
	matrix *R = new matrix();
	//init T according to input
	//multipy 
	multiply_top_stack(R);
	delete(R);
}

void rotate_x(double theta){
	matrix *R = new matrix();
	//init T according to input
	//multipy 
	multiply_top_stack(R);
	delete(R);
}

void scale(double x, double y, double z){
	matrix *S = new matrix();
	//init T according to input
	//multipy 
	multiply_top_stack(S);
	delete(S);
}

int main(){
	return 0;
}

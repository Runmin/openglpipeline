struct point3d{
	double ele[4];
};

struct point2d{
	double ele[2];
};

struct matrix{
	double ele[4][4];
};

struct linequation{
	double a,b,c;
	double x1,x2;//range of x
    double y1,y2;//range of y
};

struct planequation{
	double a,b,c,d; // do we need to have a range here?
};

struct z_buffer{
	double **ele; // initialized to image place size
};

struct stack{
	struct matrix *mptr;
	int head;
};


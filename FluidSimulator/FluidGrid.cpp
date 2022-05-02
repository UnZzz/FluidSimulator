#include "FluidGrid.h"
#include <math.h>

float visc = 0.0f;
float diff = 0.0f;
float dt = 0.4f; 
float densityAdd = 50.0f;
float forceX = 0.0f;
float forceY = -7.0f;
float circleRadius = 1.0f;


inline int GetSize(int pixel_size, int length) {
	if (length % pixel_size) {
		return length / pixel_size + 1;
	}
	else return length / pixel_size;
}

inline void set_bnd(int size_x, int size_y, int b, float** x)
{
	for (int i = 1; i <= size_x-1; i++) {
		x[0][i] = b == 1 ? -x[1][i] : x[1][i];
		x[size_x + 1][i] = b == 1 ? -x[size_x][i] : x[size_x][i];

	}
	for (int j = 1; j <= size_y - 2; j++) {
		x[j][0] = b == 2 ? -x[j][1] : x[j][1];
		x[j][size_y + 1] = b == 2 ? -x[j][size_y] : x[j][size_y];
	}
	x[0][0] = 0.5 * (x[1][0] + x[0][1]);
	x[0][size_y + 1] = 0.5 * (x[1][size_y + 1] + x[0][size_y]);
	x[size_x + 1][0] = 0.5 * (x[size_x][0] + x[size_x + 1][1]);
	x[size_x + 1][size_y + 1] = 0.5 * (x[size_x][size_y + 1] + x[size_x + 1][size_y]);
}


inline void add_source(int size_x, int size_y, float** x, float** s , float dt) {
	for (int i = 0; i < size_x; i++) {
		for (int j = 0; j < size_y; j++) {
			x[i][j] += s[i][j]*dt;
		}
	}
}

inline void diffuse(int size_x, int size_y, int b , float** x, float** x0, float diff, float dt) {
	float a = dt * diff * size_x * size_y;
	for (int k = 0; k < 20; k++) {
		for (int i = 1; i <= size_x - 2; i++) {
			for (int j = 1; j <= size_y - 2; j++) {
				x[i][j] = (x0[i][j] + a * (x0[i - 1][j] + x0[i + 1][j] + x0[i][j - 1] + x0[i][j + 1]))/(1+4*a);
			
			}
		}
		set_bnd(size_x, size_y, b, x);
	}
}

inline void advect(int size_x, int size_y, int b, float** d, float** d0, float** u, float** v, float dt) {
	int i = 0, j = 0, i0 = 0, j0 = 0, i1 = 0, j1 = 0;
	float x = 0, y = 0, s0 = 0, t0 = 0, s1 = 0, t1 = 0, dt0 = 0;
	dt0 = dt * sqrt(size_x*size_y);
	for (i = 1; i <= size_x-2; i++) {
		for (j = 1; j <= size_y-2; j++) {
			x = i - dt0 * u[i][j]; 
			y = j - dt0 * v[i][j];
			if (x < 0.5) x = 0.5; 
			if (x > size_x - 1 + 0.5) 
				x = size_x - 1 + 0.5; 
			i0 = (int)x; 
			i1 = i0 + 1;
			if (y < 0.5) y = 0.5; 
			if (y > size_y - 1 + 0.5) 
				y = size_y - 1 + 0.5; 
			j0 = (int)y; 
			j1 = j0 + 1;
			s1 = x - i0; 
			s0 = 1 - s1; 
			t1 = y - j0; 
			t0 = 1 - t1;
			d[i][j] = s0 * (t0 * d0[i0][j0] + t1 * d0[i0][j1]) +
				s1 * (t0 * d0[i1][j0] + t1 * d0[i1][j1]);
		}
	}
	set_bnd(size_x, size_y, b, d);
}

inline void dens_step(int size_x, int size_y, float** x, float** x0, float** u, float** v, float diff,
	float dt)
{
	add_source(size_x,size_y, x, x0, dt);
	std::swap(x0, x); diffuse(size_x,size_y, 0, x, x0, diff, dt);
	std::swap(x0, x); advect(size_x,size_y, 0, x, x0, u, v, dt);
}

inline void project(int size_x, int size_y, float** u, float** v, float** p, float** div)
{
	int i, j, k;
	float h;
	h = 1.0 / sqrt(size_x*size_y);
	for (i = 1; i <= size_x; i++) {
		for (j = 1; j <= size_y; j++) {
			div[i][j] = -0.5 * h * (u[i + 1][j] - u[i - 1][j] +
				v[i][j + 1] - v[i][j - 1]);
			p[i][j] = 0;
		}
	}
	set_bnd(size_x,size_y, 0, div); set_bnd(size_x,size_y, 0, p);
	for (k = 0; k < 20; k++) {
		for (i = 1; i <= size_x-2; i++) {
			for (j = 1; j <= size_y-2; j++) {
				p[i][j] = (div[i][j] + p[i - 1][ j] + p[i + 1][j] +
					p[i][j - 1] + p[i][j + 1]) / 4;
			}
		}
		set_bnd(size_x,size_y, 0, p);
	}
	for (i = 1; i <= size_x-2; i++) {
		for (j = 1; j <= size_y-2; j++) {
			u[i][j] -= 0.5 * (p[i + 1][j] - p[i - 1][j]) / h;
			v[i][j] -= 0.5 * (p[i][j + 1] - p[i][j - 1]) / h;
		}
	}
	set_bnd(size_x,size_y, 1, u); set_bnd(size_x,size_y, 2, v);
}

inline void vel_step(int size_x,int size_y, float** u, float** v, float** u0, float** v0,
	float visc, float dt)
{
	add_source(size_x,size_y, u, u0, dt); 
	add_source(size_x,size_y, v, v0, dt);
	std::swap(u0, u); diffuse(size_x, size_y, 1, u, u0, visc, dt);
	std::swap(v0, v); diffuse(size_x, size_y, 2, v, v0, visc, dt);
	project(size_x, size_y, u, v, u0, v0);
	std::swap(u0, u); 
	std::swap(v0, v);
	advect(size_x, size_y, 1, u, u0, u0, v0, dt); 
	advect(size_x, size_y, 2, v, v0, u0, v0, dt);
	project(size_x, size_y, u, v, u0, v0);
}


FluidGrid::FluidGrid(SDL_Window* window, int pixel_size, int width, int height)
	:renderer_(SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED)),
	Size_X_(GetSize(pixel_size, width)),
	Size_Y_(GetSize(pixel_size, height)),
	pixel_size_(pixel_size)
{
	Init(Size_X_+2,Size_Y_+2,pixel_size);
	return;
}

void FluidGrid::Init(int size_x,int size_y,int pixel_size)
{
	display_grid_ = new SDL_Rect * [size_x];
	for (int i = 0; i < size_x; i++) {
		display_grid_[i] = new SDL_Rect[size_y];
	}

	dens_ = new float* [size_x];
	for (int i = 0; i < size_x; i++) {
		dens_[i] = new float[size_y];
	}

	dens_prev_ = new float* [size_x];
	for (int i = 0; i < size_x; i++) {
		dens_prev_[i] = new float[size_y];
	}

	u_ = new float* [size_x];
	for (int i = 0; i < size_x; i++) {
		u_[i] = new float[size_y];
	}

	u_prev_ = new float* [size_x];
	for (int i = 0; i < size_x; i++) {
		u_prev_[i] = new float[size_y];
	}

	v_ = new float* [size_x];
	for (int i = 0; i < size_x; i++) {
		v_[i] = new float[size_y];
	}

	v_prev_ = new float* [size_x];
	for (int i = 0; i < size_x; i++) {
		v_prev_[i] = new float[size_y];
	}

	s_ = new float* [size_x];
	for (int i = 0; i < size_x; i++) {
		s_[i] = new float[size_y];
	}


	for (int i = 0; i < size_x; i++) {
		for (int j = 0; j < size_y; j++) {
			display_grid_[i][j].x = i * pixel_size;
			display_grid_[i][j].y = j * pixel_size;
			display_grid_[i][j].w = pixel_size;
			display_grid_[i][j].h = pixel_size;
			dens_[i][j] = 0;
			dens_prev_[i][j] = 0;
			v_[i][j] = 0;
			v_prev_[i][j] = 0;
			u_[i][j] = 0;
			u_prev_[i][j] = 0;
		}
	}

}

FluidGrid::~FluidGrid() {

}


void FluidGrid::Put(int x, int y) {
	int posX, posY;
	if (x % pixel_size_) {
		posX = x / pixel_size_ + 1;
	}
	else {
		posX = x / pixel_size_;
	}
	if (y % pixel_size_) {
		posY = y / pixel_size_ + 1;
	}
	else {
		posY = y / pixel_size_;
	}
	dens_prev_[posX][posY] += 0.1;
	u_prev_[posX][posY] += forceX * dt;
	v_prev_[posX][posY] += forceY * dt;
}

void FluidGrid::Quit() {
	delete display_grid_;
	delete dens_;
	delete dens_prev_;
	delete u_;
	delete u_prev_;
	delete v_;
	delete v_prev_;
	SDL_DestroyRenderer(renderer_);
}

void FluidGrid::Update() {
	vel_step(Size_X_, Size_Y_, u_, v_, u_prev_, v_prev_, visc, dt);
	dens_step(Size_X_, Size_Y_, dens_, dens_prev_, u_, v_, diff, dt);
}

void FluidGrid::Display() {
	SDL_SetRenderDrawColor(renderer_, 0, 0, 0, 0);

	SDL_RenderClear(renderer_);

	for (int i = 0; i < Size_X_; i++) {
		for (int j = 0; j < Size_Y_; j++) {
			SDL_SetRenderDrawColor(renderer_, (int)dens_[i][j] % 255, (int)dens_[i][j]%255, (int)dens_[i][j] % 255, 255);
			SDL_RenderFillRect(renderer_, &display_grid_[i][j]);
		}
	}

	SDL_RenderPresent(renderer_);

}
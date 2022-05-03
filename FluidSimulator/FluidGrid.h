#pragma once
#include <SDL.h>
#include <vector>
#include <memory>
class FluidGrid
{
public:
	FluidGrid() = default;
	FluidGrid(SDL_Window* window,int pixel_size, int N);
	~FluidGrid();
	void Update();
	void Display();
	void Quit();
	void Put(int x, int y);
	void Draw(int x, int y);
private:
	bool put;
	int x, y;

	int IX(int i, int j);
	void addSource(float* x, float* s, float dt);
	void diffuse(int N, int b, float* x, float* x0, float diff, float dt);
	void advect(int N, int b, float* d, float* d0, float* u, float* v, float dt);
	void densityStep(int N, float* x, float* x0, float* u, float* v, float diff, float dt);
	void project(int N, float* u, float* v, float* p, float* div);
	void velocityStep(int N, float* u, float* v, float* u0, float* v0, float visc, float dt);
	void setBounds(int N, int b, float* x);

	int N_;
	int pixel_size_;
	SDL_Rect* display_grid_;
	float* u_;
	float* u_prev_;
	float* v_;
	float* v_prev_;
	float* dens_;
	float* dens_prev_;
	float* s_;
	SDL_Renderer* renderer_;
};


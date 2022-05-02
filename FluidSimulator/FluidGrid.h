#pragma once
#include <SDL.h>
#include <vector>
#include <memory>
class FluidGrid
{
public:
	FluidGrid() = default;
	FluidGrid(SDL_Window* window,int pixel_size, int width, int height);
	void Init(int size_x, int size_y, int pixel_size);
	~FluidGrid();
	void Update();
	void Display();
	void Quit();
	void Put(int x, int y);
private:
	int Size_X_, Size_Y_;
	int pixel_size_;
	SDL_Rect** display_grid_;
	float** u_;
	float** u_prev_;
	float** v_;
	float** v_prev_;
	float** dens_;
	float** dens_prev_;
	float** s_;
	SDL_Renderer* renderer_;
};


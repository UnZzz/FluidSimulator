#pragma once
#include "FluidGrid.h"
#include <SDL.h>
#include <vector>
#include <memory>
class Simulator
{
public:
	Simulator(int pixel_size,int N);
	~Simulator();
	void Init();
	void Event();
	void Update();
	void Display();
	void Quit();
	bool quit;
private:
	FluidGrid fluid_grid_;
	SDL_Window* window_ ;
};


#include "Simulator.h"
#include <iostream>
Simulator::Simulator(int pixel_size, int N):
	quit(false)
{	
	if (SDL_Init(SDL_INIT_VIDEO)) {
		SDL_Log("SDL¿â³õÊ¼»¯Ê§°Ü: %s", SDL_GetError());
		return;
	}
	window_ = SDL_CreateWindow("FluidSimulator", 100, 100, N, N, 0);
	if (window_ == nullptr) {
		SDL_Log("´°¿Ú³õÊ¼»¯Ê§°Ü: %s", SDL_GetError());
		return;
	}
	fluid_grid_ = FluidGrid(window_, pixel_size, N);
}

Simulator::~Simulator() {
	SDL_DestroyWindow(window_);
	SDL_Quit();
}

void Simulator::Init() {

}

void Simulator::Event() {
	SDL_Event event;
	while (SDL_PollEvent(&event))
	{
		if (event.type == SDL_QUIT) {
			Quit();
			break;
		}


		static bool draw = false;
		if (event.type == SDL_MOUSEBUTTONDOWN) {
			draw = true;
		}
		if (event.type == SDL_MOUSEBUTTONUP) {
			draw = false;
		}
		if (draw) {
			int x = event.button.x;
			int y = event.button.y;
			fluid_grid_.Put(x, y);
		}
	}
}

void Simulator::Update() {
	fluid_grid_.Update();
}

void Simulator::Display() {
	fluid_grid_.Display();
}

void Simulator::Quit() {
	fluid_grid_.Quit();
	quit = true;
	SDL_DestroyWindow(window_);
	SDL_Quit();
}

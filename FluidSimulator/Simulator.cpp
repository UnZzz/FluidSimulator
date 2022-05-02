#include "Simulator.h"
#include <iostream>
Simulator::Simulator(int pixel_size, int width, int height):
	quit(false)
{	
	if (SDL_Init(SDL_INIT_VIDEO)) {
		SDL_Log("SDL库初始化失败: %s", SDL_GetError());
		return;
	}
	window_ = SDL_CreateWindow("FluidSimulator", 400, 200, width, height, 0);
	if (window_ == nullptr) {
		SDL_Log("窗口初始化失败: %s", SDL_GetError());
		return;
	}

	fluid_grid_ = FluidGrid(window_, pixel_size, width, height);
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
		if (event.type == SDL_MOUSEBUTTONDOWN) {
			int x = event.button.x;
			int y = event.button.y;
			std::cout << "x:" << x << ' ' << "y:" << y << std::endl;
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

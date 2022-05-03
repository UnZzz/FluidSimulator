#include "FluidGrid.h"
#include <iostream>
#include <math.h>

float visc = 0.0f;
float diff = 0.0f;
float dt = 0.5f; 
float densityAdd = 50.0f;
float forceX = 0.0f;
float forceY = -0.5f;


inline int size(int N) {
    return N * N;
}

inline int GetSize(int pixel_size, int length) {
	if (length % pixel_size) {
		return length / pixel_size + 1;
	}
	else return length / pixel_size;
}

int FluidGrid::IX(int i, int j)
{
    return i + (N_ + 2) * j;
}

void FluidGrid::addSource(float* x, float* s, float dt)
{
    for (int index = 0; index < (N_ + 2) * (N_ + 2); index++)
    {
        x[index] += dt * s[index];
    }
}

void FluidGrid::diffuse(int N, int b, float* x, float* x0, float diff, float dt)
{
    int i, j, k;
    float a = dt * diff * N * N;
    for (k = 0; k < 20; k++)
    {
        for (i = 1; i <= N; i++)
        {
            for (j = 1; j <= N; j++)
            {
                x[IX(i, j)] = (x0[IX(i, j)] + a * (x[IX(i - 1, j)] + x[IX(i + 1, j)] + x[IX(i, j - 1)] + x[IX(i, j + 1)])) / (1 + 4 * a);
            }
        }
        setBounds(N, b, x);
    }
}

void FluidGrid::advect(int N, int b, float* d, float* d0, float* u, float* v, float dt)
{
    int i, j, i0, j0, i1, j1;
    float x, y, s0, t0, s1, t1;
    float dt0 = dt * N;

    for (i = 1; i <= N; i++)
    {
        for (j = 1; j <= N; j++)
        {

            x = i - dt0 * u[IX(i, j)];
            y = j - dt0 * v[IX(i, j)];

            if (x < 0.5)
            {
                x = 0.5;
            }
            if (x > N + 0.5)
            {
                x = N + 0.5;
            }

            i0 = (int)x;
            i1 = i0 + 1;

            if (y < 0.5)
            {
                y = 0.5;
            }
            if (y > N + 0.5)
            {
                y = N + 0.5;
            }

            j0 = (int)y;
            j1 = j0 + 1;


            s1 = x - i0;
            s0 = 1 - s1;
            t1 = y - j0;
            t0 = 1 - t1;


            d[IX(i, j)] = s0 * (t0 * d0[IX(i0, j0)] + t1 * d0[IX(i0, j1)]) + s1 * (t0 * d0[IX(i1, j0)] + t1 * d0[IX(i1, j1)]);
        }
    }

    setBounds(N, b, d);
}

void FluidGrid::densityStep(int N, float* x, float* x0, float* u, float* v, float diff, float dt)
{
    float* tmp;
    addSource(x, x0, dt);
    tmp = x0;
    x0 = x;
    x = tmp;
    diffuse(N, 0, x, x0, diff, dt);
    tmp = x0;
    x0 = x;
    x = tmp; 
    advect(N, 0, x, x0, u, v, dt);
}

void FluidGrid::project(int N, float* u, float* v, float* p, float* div)
{
    int i, j, k;
    float h = 1.0f / N;

    for (i = 1; i <= N; i++)
    {
        for (j = 1; j <= N; j++)
        {
            div[IX(i, j)] = -0.5 * h * (u[IX(i + 1, j)] - u[IX(i - 1, j)] + v[IX(i, j + 1)] - v[IX(i, j - 1)]);
            p[IX(i, j)] = 0;
        }
    }

    setBounds(N, 0, div);
    setBounds(N, 0, p);

    for (k = 0; k < 20; k++)
    {
        for (i = 1; i <= N; i++)
        {
            for (j = 1; j <= N; j++)
            {
                p[IX(i, j)] = (div[IX(i, j)] + p[IX(i - 1, j)] + p[IX(i + 1, j)] + p[IX(i, j - 1)] + p[IX(i, j + 1)]) / 4.0f;
            }
        }
        setBounds(N, 0, p);
    }

    for (i = 1; i <= N; i++)
    {
        for (j = 1; j <= N; j++)
        {
            u[IX(i, j)] -= 0.5 * (p[IX(i + 1, j)] - p[IX(i - 1, j)]) / h;
            v[IX(i, j)] -= 0.5 * (p[IX(i, j + 1)] - p[IX(i, j - 1)]) / h;
        }
    }
    setBounds(N, 1, u);
    setBounds(N, 2, v);
}

void FluidGrid::velocityStep(int N, float* u, float* v, float* u0, float* v0, float visc, float dt)
{
    float* tmp;

    addSource(u, u0, dt);
    addSource(v, v0, dt);

    tmp = u0;
    u0 = u;
    u = tmp; 
    diffuse(N, 1, u, u0, visc, dt);

    tmp = v0;
    v0 = v;
    v = tmp; 
    diffuse(N, 2, v, v0, visc, dt);
    project(N, u, v, u0, v0);

    tmp = u0;
    u0 = u;
    u = tmp; 
    tmp = v0;
    v0 = v;
    v = tmp; 
    advect(N, 1, u, u0, u0, v0, dt);
    advect(N, 2, v, v0, u0, v0, dt);
    project(N, u, v, u0, v0);
}

void FluidGrid::setBounds(int N, int b, float* x)
{

    for (int i = 1; i <= N; i++)
    {
        x[IX(0, i)] = b == 1 ? -x[IX(1, i)] : x[IX(1, i)];
        x[IX(N + 1, i)] = b == 1 ? -x[IX(N, i)] : x[IX(N, i)];
        x[IX(i, 0)] = b == 2 ? -x[IX(i, 1)] : x[IX(i, 1)];
        x[IX(i, N + 1)] = b == 2 ? -x[IX(i, N)] : x[IX(i, N)];
    }


    x[IX(0, 0)] = 0.5 * (x[IX(1, 0)] + x[IX(0, 1)]);
    x[IX(0, N + 1)] = 0.5 * (x[IX(1, N + 1)] + x[IX(0, N)]);
    x[IX(N + 1, 0)] = 0.5 * (x[IX(N, 0)] + x[IX(N + 1, 1)]);
    x[IX(N + 1, N + 1)] = 0.5 * (x[IX(N, N + 1)] + x[IX(N + 1, N)]);
}


FluidGrid::FluidGrid(SDL_Window* window, int pixel_size, int N)
	:renderer_(SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED)),
	N_(GetSize(pixel_size, N)),
	pixel_size_(pixel_size),
    put(false)
{	
    display_grid_ = new SDL_Rect[size(N_ + 2)];
    u_ = new float[size(N_ + 2)];
	v_ = new float[size(N_ + 2)];
	u_prev_ = new float[size(N_ + 2)];
	v_prev_ = new float[size(N_ + 2)];
	dens_ = new float[size(N_ + 2)];
	dens_prev_ = new float[size(N_ + 2)];

    for (int i = 0; i < N_ + 2; i++) {
        for (int j = 0; j < N_ + 2; j++) {

            u_[IX(i, j)] = 0;
            v_[IX(i, j)] = 0;
            
            u_prev_[IX(i, j)] = 0;
            v_prev_[IX(i, j)] = 0;

            dens_[IX(i, j)] = 0;
            dens_prev_[IX(i, j)] = 0;

            display_grid_[IX(i, j)].x = i * pixel_size_;
            display_grid_[IX(i, j)].y = j * pixel_size_;
            display_grid_[IX(i, j)].w = pixel_size_;
            display_grid_[IX(i, j)].h = pixel_size_;
        }
    }

	return;
}


FluidGrid::~FluidGrid() {

}

void FluidGrid::Draw(int x, int y) {
    static int x_prev = 0, y_prev = 0;
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

    int q = 2;
    
    for (int i = posX - q; i <= posX + q; i++) {
        for (int j = posY - q; j <= posY + q; j++) {
            if (!(sqrt((x - x_prev) * (x - x_prev) + (y - y_prev) * (y - y_prev)) > 100)) {
                dens_prev_[IX(i, j)] += 200;
                u_prev_[IX(i, j)] += (float)(x - x_prev) * dt * 0.1;
                v_prev_[IX(i, j)] += (float)(y - y_prev) * dt * 0.1;
            }
            
        }
    }
    x_prev = x;
    y_prev = y;
}

void FluidGrid::Put(int x, int y) {
    put = true;
    this->x = x;
    this->y = y;
}

void FluidGrid::Quit() {
    std::cout << "Destroying FluidGrid" << std::endl;
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
    for (int i = 0; i < size(N_ + 2); i++) {
            dens_prev_[i] = 0;
            u_prev_[i] = 0;
            v_prev_[i] = 0;
    }


    if (put) {
        std::cout << x<< " " << y <<std::endl;
        Draw(x, y);
        put = false;
    }
    velocityStep(N_, u_, v_, u_prev_, v_prev_, visc, dt);
    densityStep(N_, dens_, dens_prev_, u_, v_, diff, dt);
}

inline float get_R(float H, float S, float V) {
    float R;
    if (S == 0) {
        R = V;
    }
    else {
        H /= 60;
        int i = H;
        float f = H - i;
        float a = V * (1 - S);
        float b = V * (1 - S * f);
        float c = V * (1 - S * (1 - f));
        if (i == 0) {
            R = V;
        }
        if (i == 1) {
            R = b;
        }
        if (i == 2) {
            R = a;
        }
        if (i == 3) {
            R = a;
        }
        if (i == 4) {
            R = c;
        }
        if (i == 5) {
            R = V;
        }
    }
    return R*255;
}

inline float get_G(float H, float S, float V) {
    float G;
    if (S == 0) {
        G = V;
    }
    else {
        H /= 60;
        int i = H;
        float f = H - i;
        float a = V * (1 - S);
        float b = V * (1 - S * f);
        float c = V * (1 - S * (1 - f));
        if (i == 0) {
            G = c;
        }
        if (i == 1) {
            G = V;
        }
        if (i == 2) {
            G = V;
        }
        if (i == 3) {
            G = b;
        }
        if (i == 4) {
            G = a;
        }
        if (i == 5) {
            G = a;
        }
    }
    return G*255;
}

inline float get_B(float H, float S, float V) {
    float B;
    if (S == 0) {
        B = V;
    }
    else {
        H /= 60;
        int i = H;
        float f = H - i;
        float a = V * (1 - S);
        float b = V * (1 - S * f);
        float c = V * (1 - S * (1 - f));
        if (i == 0) {
            B = a;
        }
        if (i == 1) {
            B = a;
        }
        if (i == 2) {
            B = c;
        }
        if (i == 3) {
            B = V;
        }
        if (i == 4) {
            B = V;
        }
        if (i == 5) {
            B = b;
        }
    }
    return B*255;
}

void FluidGrid::Display() {

    SDL_SetRenderDrawColor(renderer_, 0, 0, 0, 0);

	SDL_RenderClear(renderer_);

    for (int i = 0; i < size(N_ + 2); i++) {
        SDL_SetRenderDrawColor(renderer_,
            get_B((int)dens_[i] % 360, 1, dens_[i] > 25 ? 1 : 0),
            get_R((int)dens_[i] % 360, 0.5, dens_[i] > 25 ? 1 : 0),
            get_G((int)dens_[i] % 360, 1, dens_[i] > 25 ? 1 : 0),
            255);
        SDL_RenderFillRect(renderer_, display_grid_+i);
    }

	SDL_RenderPresent(renderer_);
}
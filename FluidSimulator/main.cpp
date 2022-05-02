#include "Simulator.h"
int main(int argc, char* argv[]) {
	Simulator simulator(4,640,480);
	while (true) {
		simulator.Event();
		if (simulator.quit)
			break;
		simulator.Update();
		simulator.Display();
	}
	return 0;
}
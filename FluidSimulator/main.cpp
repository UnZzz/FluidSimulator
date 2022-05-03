#include "Simulator.h"
int main(int argc, char* argv[]) {
	Simulator simulator(8,640);
	while (true) {
		simulator.Event();
		if (simulator.quit)
			break;
		simulator.Update();
		simulator.Display();
	}
	return 0;
}
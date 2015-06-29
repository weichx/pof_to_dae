#pragma once
#include <string>
#include <iostream>
using namespace std;

class AsyncProgress {

public:

	AsyncProgress() : target(0), progress(0), interrupted(false), terminated(false), errorcode(0), finished(false) {}

	unsigned int target;
	unsigned int progress;
	bool interrupted;
	bool terminated;
	std::string message;
	int errorcode;
	bool finished;
	bool doYieldOnNotify;

	void setTarget(unsigned int target) { this->target = target; }
	unsigned int GetTarget() { return this->target; }

	void incrementWithMessage(std::string msg) {
		this->progress++; 
		this->setMessage(msg);
		this->Notify(); 
	}
	void incrementProgress() { 
		this->progress++;
		this->Notify();
	}

	void Notify() {
		cout << this->message.c_str() << endl;
	}

	float GetTargetPercent() { return 100.00f * (float(this->progress) / float(this->target)); }

	void SetProgress(unsigned int progress) { this->progress = progress; }
	unsigned int GetProgress() { return this->progress; }
	void IncrementProgress() { this->progress++; }

	void Interrupt() { this->interrupted = true; }
	bool Interrupted() { return this->interrupted; }

	void EarlyTerminate() { this->terminated = true; }
	bool EarlyTerminated() { return this->terminated; }

	bool Finished() { return this->finished; }
	void Finish() { this->finished = true; }

	void setMessage(char * pMessage) {
		this->message = string(pMessage);
	}

	void setMessage(std::string &msg) { this->message = msg; }
	std::string GetMessage() { return this->message; }

	void SetError(int errorcode) { this->errorcode = errorcode; }
	int GetError() { return this->errorcode; }
	virtual ~AsyncProgress() { }
};
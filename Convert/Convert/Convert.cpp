// Convert.cpp : Defines the entry point for the console application.
//
#include <tchar.h>
#include <iostream>
#include "Src/pcs_file.h"


void SaveModel(std::string path, PCS_Model * pModel, AsyncProgress * pAsyncProgress) {
	int err = pModel->SaveToDAE(path, pAsyncProgress, 1, 0);
	if (err != 0)
	{
		switch (err)
		{
		case 1:
			pAsyncProgress->setMessage("Collada Save Failed");
			break;

		default:
			pAsyncProgress->setMessage("No error message exists for this error, bitch at kazan!");
			break;
		}
	}
}

int _tmain(int argc, _TCHAR* argv[])
{
	PCS_Model * pModel = new PCS_Model();
	AsyncProgress *pAsyncProgress = new AsyncProgress();
	char name[256];
	pModel->LoadFromPOF("C:\\Users\\Matt\\Desktop\\Freespace Models\\Damocles\\Damocles.pof", pAsyncProgress);
	SaveModel("C:\\Users\\Matt\\Desktop\\Freespace Models\\Damocles\\Damocles.xml", pModel, pAsyncProgress);
	delete pModel;
	delete pAsyncProgress;
	std::cin.getline(name, 256);
	return 0;
}


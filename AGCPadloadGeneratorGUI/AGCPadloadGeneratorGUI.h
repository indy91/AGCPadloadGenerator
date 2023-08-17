
// AGCPadloadGeneratorGUI.h:
//

#pragma once

#ifndef __AFXWIN_H__
	#error "include 'pch.h' before including this file for PCH"
#endif

#include "resource.h"


// CAGCPadloadGeneratorGUIApp
//

class CAGCPadloadGeneratorGUIApp : public CWinApp
{
public:
	CAGCPadloadGeneratorGUIApp();

public:
	virtual BOOL InitInstance();

	DECLARE_MESSAGE_MAP()
};

//extern CAGCPadloadGeneratorGUIApp theApp;

namespace Utilities
{
	double Text2Double(CEdit *ed);
	int Text2Octal(CEdit *ed);
	void Double2Text(double val, CEdit *ed, int length);
}
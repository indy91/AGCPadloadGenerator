// MainMenuDialog.cpp
//

#include "pch.h"
#include "AGCPadloadGeneratorGUI.h"
#include "MainMenuDialog.h"
#include "afxdialogex.h"
#include "AGCPadloadGeneratorGUIDlg.h"
#include "LGCPadloadGenerator.h"
#include "BlockIPadloadGenerator.h"
#include "SkylarkPadloadGenerator.h"

// MainMenuDialog-Dialog

IMPLEMENT_DYNAMIC(MainMenuDialog, CDialogEx)

MainMenuDialog::MainMenuDialog(CWnd* pParent /*=nullptr*/)
	: CDialogEx(IDD_DIALOG1, pParent)
{
	cmc_form = NULL;
	lgc_form = NULL;
	blockI_form = NULL;
	skylark_form = NULL;
}

MainMenuDialog::~MainMenuDialog()
{
	if (cmc_form)
	{
		delete cmc_form;
		cmc_form = NULL;
	}
	if (lgc_form)
	{
		delete lgc_form;
		lgc_form = NULL;
	}
	if (blockI_form)
	{
		delete blockI_form;
		blockI_form = NULL;
	}
	if (skylark_form)
	{
		delete skylark_form;
		skylark_form = NULL;
	}
}

void MainMenuDialog::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
}


BEGIN_MESSAGE_MAP(MainMenuDialog, CDialogEx)
	ON_BN_CLICKED(IDC_BUTTON1, &MainMenuDialog::OnBnClickedCMC)
	ON_BN_CLICKED(IDC_BUTTON2, &MainMenuDialog::OnBnClickedLGC)
	ON_BN_CLICKED(IDC_BUTTON3, &MainMenuDialog::OnBnClickedBlockI)
	ON_BN_CLICKED(IDC_BUTTON4, &MainMenuDialog::OnBnClickedSkylark)
END_MESSAGE_MAP()


void MainMenuDialog::OnBnClickedCMC()
{
	if (cmc_form == NULL)
	{
		cmc_form = new CAGCPadloadGeneratorGUIDlg(this);
		cmc_form->DoModal();
	}
	//cmc_form->ShowWindow(SW_SHOW);
}


void MainMenuDialog::OnBnClickedLGC()
{
	if (lgc_form == NULL)
	{
		lgc_form = new LGCPadloadGenerator(this);
		lgc_form->DoModal();
	}
	//lgc_form->ShowWindow(SW_SHOW);
}


void MainMenuDialog::OnBnClickedBlockI()
{
	if (blockI_form == NULL)
	{
		blockI_form = new BlockIPadloadGenerator(this);
		blockI_form->DoModal();
	}
	//lgc_form->ShowWindow(SW_SHOW);
}


void MainMenuDialog::OnBnClickedSkylark()
{
	if (skylark_form == NULL)
	{
		skylark_form = new SkylarkPadloadGenerator(this);
		skylark_form->DoModal();
	}
}

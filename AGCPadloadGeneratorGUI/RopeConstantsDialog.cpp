// RopeConstantsDialog.cpp : implementation file
//

#include "pch.h"
#include "AGCPadloadGeneratorGUI.h"
#include "RopeConstantsDialog.h"
#include "afxdialogex.h"
#include "AGCPadloadGenerator.h"

// RopeConstantsDialog dialog

IMPLEMENT_DYNAMIC(RopeConstantsDialog, CDialogEx)

RopeConstantsDialog::RopeConstantsDialog(CWnd* pParent /*=nullptr*/)
	: CDialogEx(IDD_DIALOG5, pParent)
{

}

RopeConstantsDialog::~RopeConstantsDialog()
{
}

void RopeConstantsDialog::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
	DDX_Control(pDX, IDC_EDIT1, RopeConstantsYear);
	DDX_Control(pDX, IDC_EDIT2, OutputBox);
}


BEGIN_MESSAGE_MAP(RopeConstantsDialog, CDialogEx)
	ON_BN_CLICKED(IDOK, &RopeConstantsDialog::OnBnClickedGenerate)
END_MESSAGE_MAP()


// RopeConstantsDialog message handlers


void RopeConstantsDialog::OnBnClickedGenerate()
{
	AGCPadloadGenerator agc;

	int year = Utilities::Text2Int(&RopeConstantsYear);

	agc.GenerateRopeConstants(year);

	OutputBox.SetWindowTextW(L"Constants written to RopeConstants.txt!");
}

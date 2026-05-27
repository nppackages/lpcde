********************************************************************************
* LPCDE STATA PACKAGE -- build Mata library
********************************************************************************
*!version 1.0.0 2026-05-27

version 16.0
clear all
set more off

local src "stata/lpcde_fun.do"
local outdir "stata"
capture confirm file "`src'"
if (_rc) {
	local src "lpcde_fun.do"
	local outdir "."
	capture confirm file "`src'"
	if (_rc) {
		display as error "lpcde_fun.do not found"
		exit 601
	}
}

mata: mata clear
quietly do "`src'"
lmbuild lpcde, dir("`outdir'") replace
quietly adopath ++ "`outdir'"
mata: mata mlib index

display as text "Built `outdir'/lpcde.mlib"

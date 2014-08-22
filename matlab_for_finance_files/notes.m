%% Matlab notes
% This script contains notes, shortcuts and customizations for Matlab.

%% Course requirements:
% - git and github account
% - MATLAB coursework account

%% Customizations

% preferences: (leave open at beginning)
% - create own shortcuts file
% - adapt font size
% - Editor display: 
% - - highlight current line
% - - show line numbers
% - - enable datatips

% - get rid of toolbars?

% beep off?

%% MATLAB IDE features
% - Debugger
% - profiler
% - command window shortcuts: up, down, autocomplete
% - file diff
% - change variable names for multiple occurrences

%% Some additional language features
% - anonymous functions
% - x(:)
% - varfun: apply function to table columns

%% Further interesting topics
% - Quandl
% - parallel computing
% - Dataset / Table / Time Series objects
% - Unit testing

%% General programming topics
% - git
% - best practices
% - folder structure
% - floating point arithmetic

%% Shortcuts

% open new script file:
% C-x n

% next tab
% C-x b

% previous tab
% C-x C-b

% next cell
% C-M-n

% previous cell
% C-M-p

% store bookmark
% C-F2

% next bookmark
% F2


%% Notes

% folder structure: 
% - one folder per project
% - auxiliary files in folder lib
% - saved pictures in folder pics
% - saved data in folder data

% set current directory
% -> current directory is default path to save pictures, data, ...
% -> always be aware of the current directory you are in!!
% -> at beginning of session, set current directory to current project
%    folder: use Current Folder utility
% -> common utilities of several projects should be stored in one folder,
%    and included as symbolic links (relative paths!)
% -> always store your files with relative paths!!

% display current directory
cd

% set matlab path
% right-click on project folder
% -> Add to path
% -> Selected folders and subfolders


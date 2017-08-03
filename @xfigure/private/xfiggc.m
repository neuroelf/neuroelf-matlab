function xfiggc
%XFIGURE::XFIGGC  Garbage collection of XFIGURE class.

% global reference storage
global xfigmlup xfigsngl xfigures;

% remove invalid handles (deep check!)
ivh = ~isvalidxfigure(xfigures);
xfigmlup(ivh) = [];
xfigures(ivh) = [];

% reset garbage collection counter
xfigsngl.gccount = 0;

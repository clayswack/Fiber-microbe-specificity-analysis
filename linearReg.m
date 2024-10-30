function [mdl, xModel, yModel, pValTxt, R2txt] = linearReg(x, y)

    % Linear regression
    mdl = fitlm(x,y);
    xModel = linspace(min(x), max(x), 100);
    yModel = mdl.Coefficients.Estimate(2)*xModel + mdl.Coefficients.Estimate(1);
    
    % P val text
    pValTxt =  num2str(round(mdl.Coefficients.pValue(2),4));
    pValTxt = strcat(', p =', '{ }', pValTxt);
    
    % R2 text
    R2txt =  num2str(round(mdl.Rsquared.Ordinary,3));
    R2txt = strcat(' R^2 =', '{ }', R2txt);

end
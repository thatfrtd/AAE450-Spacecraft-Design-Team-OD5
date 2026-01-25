function [year] = sec_to_year(sec)
%SEC_TO_YEAR Summary of this function goes here
%   Detailed explanation goes here
year = sec / (60 * 60 * 24 * 365.25);
end
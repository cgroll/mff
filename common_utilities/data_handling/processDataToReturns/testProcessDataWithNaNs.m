function tests = testProcessDataWithNaNs
tests = functiontests(localfunctions);

end


function testCreateValidName(testCase)

% Test specific code
invalidName = '^DAXI';
validName = createValidName(invalidName);
verifyTrue(testCase, isvarname(validName))

end

function testCreateValidNameDax(testCase)

% Test specific code
invalidName = '^DAXI';
expOutName = 'DAXI';
validName = createValidName(invalidName);
verifyEqual(testCase, expOutName, validName)

end
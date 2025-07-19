
#include <real2DFFT.H>
#include <complexFFT.H>

#include <random>
#include <cmath>
#include <numbers>

#include <gtest/gtest.h>

TEST ( MFFM_FFTW, Test2D )
{
	int x=9, y=9;
	real2DFFTData *fftData = new real2DFFTData(x,y);
	real2DFFT *fft= new real2DFFT(fftData);

	// clear the data
	fftData->clearInput();
	fftData->clearOutput();

	int temp=x/2, temp2=y/2;
	for (int j=0;j<x;j++)
		fftData->in[temp2+j*x]=1.0;
	for (int j=0;j<y;j++)
		fftData->in[temp*y+j]=1.0;
	//  fftData->in[temp*y+(y-1)/2]=20000.0;

	for (int i=0;i<fftData->getXSize();i++){
		for (int j=0;j<fftData->getYSize();j++)
			cout<<fftData->in[i*x+j]<<'\t';
		cout<<endl;
	}
	cout<<'\n'<<endl;
	fft->fwdTransform();
	fftData->reScale();
	fftData->compPowerSpec();
	fft->invTransform();

	for (int i=0;i<fftData->getXSize();i++){
		for (int j=0;j<fftData->getYSize();j++)
			cout<<fftData->in[i*x+j]<<'\t';
		cout<<endl;
	}
	cout<<'\n'<<endl;
	/*  for (int i=0;i<fftData.getXSize();i++){
		for (int j=0;j<fftData.getYHalfSize();j++)
			cout<<fftData.out[i][j].im<<'\t';
		cout<<endl;
		}*/
	for (int i=0;i<x;i++){
		for (int j=0;j<y/2+1;j++)
			cout<<fftData->power[i*(y/2+1)+j]<<'\t';
		cout<<endl;
	}
	delete fftData;
	delete fft;
}

TEST ( MFFM_FFTW, Test1D_Complex )
{
	int count = 8;
	complexFFTData fftData(count);
	complexFFT fft(&fftData);

	for (int i=0;i<count;i++){
		c_re(fftData.in[i]) = static_cast<double>(i);
		c_im(fftData.in[i]) = static_cast<double>(0);
	}

	fftw_real *temp = &c_re (fftData.in[0] );
	fftw_real *temp2 = &c_im (fftData.in[0] );
	for (int i=0; i < count; i++)
		cout << temp[i] << " + " << temp2[i] << endl;

	// forward transform :
	fft.fwdTransform();
	// inverse transform :
	fft.invTransform();

	for (int i=0; i<count; i++)
		cout << c_re(fftData.in[i]) <<' '<< c_im(fftData.in[i]) << endl;

	for (int i=0; i<count; i++)
		cout << c_re(fftData.out[i]) <<' '<< c_im(fftData.out[i]) << endl;

	// Find the power spectrum ...
	fftData.compPowerSpec();
	for (int i=0; i<count; i++)
		cout << fftData.power_spectrum[i]<<endl;
}

TEST ( MFFM_FFTW, Test1D_surrogate )
{
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(-std::numbers::pi, std::numbers::pi);

	int count = 8;
	complexFFTData fftData(count);
	complexFFT fft(&fftData);

	for (int i=0;i<count;i++){
		c_re(fftData.in[i]) = static_cast<double>(i);
		c_im(fftData.in[i]) = static_cast<double>(0);
	}

	// forward transform :
	fft.fwdTransform();

	cout << "FT original" << endl;
	for (int i=0; i<count; i++)
		cout << c_re(fftData.out[i]) <<' '<< c_im(fftData.out[i]) << endl;

	// Find the power spectrum ...
	cout << "Power spectrum original" << endl;
	fftData.compPowerSpec();
	for (int i=0; i<count; i++)
		cout << fftData.power_spectrum[i]<<endl;

	for (int i=1; i < count/2; i++)
	{
		double phase = dis(gen);
		double real_shift = cos(phase);
		double imag_shift = sin(phase);
		c_re(fftData.out[i]) = sqrt( fftData.power_spectrum[i] ) * real_shift;
		c_im(fftData.out[i]) = sqrt( fftData.power_spectrum[i] ) * imag_shift;

		c_re(fftData.out[count - i]) = c_re(fftData.out[i]);
		c_im(fftData.out[count - i]) = - c_im(fftData.out[i]);
	}

	cout << "FT surrogate" << endl;
	for (int i=0; i<count; i++)
	{
		cout << c_re(fftData.out[i]) <<' '<< c_im(fftData.out[i]) << endl;
	}

	fft.invTransform();

	for (int i=0; i<count; i++)
	{
		c_re(fftData.in[i]) /= count;
		c_im(fftData.in[i]) /= count;
	}

	cout << "surrogate" << endl;
	for (int i=0; i<count; i++)
		cout << c_re(fftData.in[i]) <<' '<< c_im(fftData.in[i]) << endl;

	fft.fwdTransform();

	cout << "Power spectrum surrogate" << endl;
	fftData.compPowerSpec();
	for (int i=0; i<count; i++)
		cout << fftData.power_spectrum[i]<<endl;
}

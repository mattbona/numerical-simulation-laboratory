#ifndef ASSET_H
#define ASSET_H

class Asset{

private:
	double asset_price;
protected:

public:
	Asset();
	~Asset();

	// Methods
	virtual void UpdateAssetPrice(double t){};
	void SetAssetPrice(double my_asset_value);
	double GetAssetPrice();
};

class GBM :public Asset{
	/*
	A class that implement the geometric brownian motion with mean=mu and std.dev=sigma
	for the method UpdateAssetValue(t).
	*/
private:
	double mu, sigma;
	double gaussian_var;
protected:

public:
	GBM(double my_mu=0, double my_sigma=1);
	~GBM();

	// Methods
	void UpdateAssetPrice(double t);
	void SetGaussianVar(double my_gaussian_var);
	double GetMu();
	double GetSigma();
};

#endif // ASSET_H

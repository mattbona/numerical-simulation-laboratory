#ifndef ASSET_H
#define ASSET_H

class Asset{

private:
	double asset_value;
protected:

public:
	Asset();
	~Asset();

	// Methods
	virtual void UpdateAssetValue(double t){};
	void SetAssetValue(double my_asset_value);
	double GetAssetValue();
};

class GBM :public Asset{

private:
	double mu, sigma;
	double gaussian_var;
protected:

public:
	GBM(double my_mu=0, double my_sigma=1);
	~GBM();

	// Methods
	void UpdateAssetValue(double t);
	void SetGaussianVar(double my_gaussian_var);
};

#endif // ASSET_H

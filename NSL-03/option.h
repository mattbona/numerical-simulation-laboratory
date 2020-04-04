#ifndef OPTION_H
#define OPTION_H

class Option{

private:

protected:
	double call_option_profit,put_option_profit;

public:
	Option(){
        call_option_profit=0;
        put_option_profit=0;
    };
	~Option(){};
	// Methods
	double GetCallOptionProfit();
    void SetCallOptionProfit(double my_call_option_profit);
	virtual void UpdateCallOptionProfit(){};
	virtual void UpdatePutOptionProfit(){};
};

class European : public Option{
private:
	double asset_price;
	double strike_price;
protected:

public:
	European();
	~European();
	// Methods
	void UpdateCallOptionProfit();
	void UpdatePutOptionProfit();
    void SetStrikePrice(double my_strike_price);
    void SetAssetPrice(double my_asset_price);
};

#endif // OPTION_H

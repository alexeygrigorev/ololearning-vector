#include <string>

class Regressor {
private:
    std::string name;
public:
    Regressor();
    virtual fit display();
    virtual std::string get_name();
};

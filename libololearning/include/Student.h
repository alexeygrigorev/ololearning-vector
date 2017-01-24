#include <string>

class Student{
private:
    std::string name;
public:
    Student(std::string);
    virtual void display();
    virtual std::string get_name();
};

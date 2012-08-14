/**
 * @file OptionSet.h
 * @brief OptionSet Class for dealing with program options
 * @author Yu Peng
 * @version 0.18
 * @date 2010-10-26
 */

#ifndef __OPTION_SET_H_

#define __OPTION_SET_H_

#include "globals.h"

#include <string>
#include <map>
#include <vector>


class OptionSet
{
public:
    enum OptionType { OptionBool, OptionInt, OptionDouble, OptionString };

    struct Option
    {
        Option(const std::string &long_name, const std::string &short_name,
                void *pointer, OptionType type, const std::string &description, 
                const std::string &default_value);

        void Process();
        std::string ToString();

        std::string long_name;
        std::string short_name;
        void *pointer;
        OptionType type;
        std::string description;
        std::string default_value;
        std::string value;
    };

	void AddOption(const std::string &long_name, const std::string &short_name,
		bool &bool_option, const std::string &description);

	void AddOption(const std::string &long_name, const std::string &short_name,
		int &int_option, const std::string &description);

	void AddOption(const std::string &long_name, const std::string &short_name,
		double &double_option, const std::string &description);

	void AddOption(const std::string &long_name, const std::string &short_name,
		std::string &string_option, const std::string &description);

	void ProcessOptions(int &argc, char *argv[]);

    std::string ToString();

private:
	void AddOption(const std::string &long_name, const std::string &short_name,
		void *pointer, OptionType type, const std::string &description, const std::string &default_value);

	std::map<std::string, void *> parameters;
	std::map<std::string, OptionType> types;
	std::vector<Option> options;
};

#endif


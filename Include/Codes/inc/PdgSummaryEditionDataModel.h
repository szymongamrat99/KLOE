//  To parse this JSON data, first install
//
//      Boost     http://www.boost.org
//      json.hpp  https://github.com/nlohmann/json
//
//  Then include this file, and then do
//
//     SummaryEdition data = nlohmann::json::parse(jsonString);

#pragma once

#include <boost/optional.hpp>
#include "nlohmann/json.hpp"

#include <boost/optional.hpp>
#include <stdexcept>
#include <regex>

#ifndef NLOHMANN_OPT_HELPER
#define NLOHMANN_OPT_HELPER
namespace nlohmann {
    template <typename T>
    struct adl_serializer<std::shared_ptr<T>> {
        static void to_json(json & j, const std::shared_ptr<T> & opt) {
            if (!opt) j = nullptr; else j = *opt;
        }

        static std::shared_ptr<T> from_json(const json & j) {
            if (j.is_null()) return std::make_shared<T>(); else return std::make_shared<T>(j.get<T>());
        }
    };
    template <typename T>
    struct adl_serializer<boost::optional<T>> {
        static void to_json(json & j, const boost::optional<T> & opt) {
            if (!opt) j = nullptr; else j = *opt;
        }

        static boost::optional<T> from_json(const json & j) {
            if (j.is_null()) return boost::optional<T>(); else return boost::optional<T>(j.get<T>());
        }
    };
}
#endif

namespace PDGProvider {
    using nlohmann::json;

    #ifndef NLOHMANN_UNTYPED_PDGProvider_HELPER
    #define NLOHMANN_UNTYPED_PDGProvider_HELPER
    inline json get_untyped(const json & j, const char * property) {
        if (j.find(property) != j.end()) {
            return j.at(property).get<json>();
        }
        return json();
    }

    inline json get_untyped(const json & j, std::string property) {
        return get_untyped(j, property.data());
    }
    #endif

    #ifndef NLOHMANN_OPTIONAL_PDGProvider_HELPER
    #define NLOHMANN_OPTIONAL_PDGProvider_HELPER
    template <typename T>
    inline std::shared_ptr<T> get_heap_optional(const json & j, const char * property) {
        auto it = j.find(property);
        if (it != j.end() && !it->is_null()) {
            return j.at(property).get<std::shared_ptr<T>>();
        }
        return std::shared_ptr<T>();
    }

    template <typename T>
    inline std::shared_ptr<T> get_heap_optional(const json & j, std::string property) {
        return get_heap_optional<T>(j, property.data());
    }
    template <typename T>
    inline boost::optional<T> get_stack_optional(const json & j, const char * property) {
        auto it = j.find(property);
        if (it != j.end() && !it->is_null()) {
            return j.at(property).get<boost::optional<T>>();
        }
        return boost::optional<T>();
    }

    template <typename T>
    inline boost::optional<T> get_stack_optional(const json & j, std::string property) {
        return get_stack_optional<T>(j, property.data());
    }
    #endif

    enum class Type : int { BEST_LIMIT, OUR_AVERAGE, OUR_FIT, HARDCODED };

    class PdgValue {
        public:
        PdgValue() = default;
        virtual ~PdgValue() = default;

        private:
        boost::optional<double> value;
        boost::optional<double> error_positive;
        boost::optional<double> error_negative;
        std::string value_text;
        Type type;
        boost::optional<double> scale_factor;
        boost::optional<bool> is_limit;
        boost::optional<bool> is_upper_limit;
        boost::optional<int64_t> confidence_level;
        boost::optional<std::string> unit;
        boost::optional<std::string> comment;

        public:
        const boost::optional<double> & get_value() const { return value; }
        boost::optional<double> & get_mutable_value() { return value; }
        void set_value(const boost::optional<double> & value) { this->value = value; }

        boost::optional<double> get_error_positive() const { return error_positive; }
        void set_error_positive(boost::optional<double> value) { this->error_positive = value; }

        boost::optional<double> get_error_negative() const { return error_negative; }
        void set_error_negative(boost::optional<double> value) { this->error_negative = value; }

        const std::string & get_value_text() const { return value_text; }
        std::string & get_mutable_value_text() { return value_text; }
        void set_value_text(const std::string & value) { this->value_text = value; }

        const Type & get_type() const { return type; }
        Type & get_mutable_type() { return type; }
        void set_type(const Type & value) { this->type = value; }

        boost::optional<double> get_scale_factor() const { return scale_factor; }
        void set_scale_factor(boost::optional<double> value) { this->scale_factor = value; }

        boost::optional<bool> get_is_limit() const { return is_limit; }
        void set_is_limit(boost::optional<bool> value) { this->is_limit = value; }

        boost::optional<bool> get_is_upper_limit() const { return is_upper_limit; }
        void set_is_upper_limit(boost::optional<bool> value) { this->is_upper_limit = value; }

        boost::optional<int64_t> get_confidence_level() const { return confidence_level; }
        void set_confidence_level(boost::optional<int64_t> value) { this->confidence_level = value; }

        boost::optional<std::string> get_unit() const { return unit; }
        void set_unit(boost::optional<std::string> value) { this->unit = value; }

        boost::optional<std::string> get_comment() const { return comment; }
        void set_comment(boost::optional<std::string> value) { this->comment = value; }
    };

    class BranchingFraction {
        public:
        BranchingFraction() = default;
        virtual ~BranchingFraction() = default;

        private:
        boost::optional<std::string> pdgid;
        boost::optional<std::string> description;
        boost::optional<int64_t> mode_number;
        boost::optional<std::string> section;
        boost::optional<std::vector<PdgValue>> pdg_values;

        public:
        boost::optional<std::string> get_pdgid() const { return pdgid; }
        void set_pdgid(boost::optional<std::string> value) { this->pdgid = value; }

        boost::optional<std::string> get_description() const { return description; }
        void set_description(boost::optional<std::string> value) { this->description = value; }

        boost::optional<int64_t> get_mode_number() const { return mode_number; }
        void set_mode_number(boost::optional<int64_t> value) { this->mode_number = value; }

        boost::optional<std::string> get_section() const { return section; }
        void set_section(boost::optional<std::string> value) { this->section = value; }

        boost::optional<std::vector<PdgValue>> get_pdg_values() const { return pdg_values; }
        void set_pdg_values(boost::optional<std::vector<PdgValue>> value) { this->pdg_values = value; }
    };

    class Property {
        public:
        Property() = default;
        virtual ~Property() = default;

        private:
        boost::optional<std::string> pdgid;
        boost::optional<std::string> description;
        boost::optional<std::vector<PdgValue>> pdg_values;

        public:
        boost::optional<std::string> get_pdgid() const { return pdgid; }
        void set_pdgid(boost::optional<std::string> value) { this->pdgid = value; }

        boost::optional<std::string> get_description() const { return description; }
        void set_description(boost::optional<std::string> value) { this->description = value; }

        boost::optional<std::vector<PdgValue>> get_pdg_values() const { return pdg_values; }
        void set_pdg_values(boost::optional<std::vector<PdgValue>> value) { this->pdg_values = value; }
    };

    class Summaries {
        public:
        Summaries() = default;
        virtual ~Summaries() = default;

        private:
        boost::optional<std::vector<Property>> properties;
        boost::optional<std::vector<BranchingFraction>> branching_fractions;

        public:
        boost::optional<std::vector<Property>> get_properties() const { return properties; }
        void set_properties(boost::optional<std::vector<Property>> value) { this->properties = value; }

        boost::optional<std::vector<BranchingFraction>> get_branching_fractions() const { return branching_fractions; }
        void set_branching_fractions(boost::optional<std::vector<BranchingFraction>> value) { this->branching_fractions = value; }
    };

    class SummaryEdition {
        public:
        SummaryEdition() = default;
        virtual ~SummaryEdition() = default;

        private:
        boost::optional<int64_t> status_code;
        boost::optional<std::string> status_message;
        boost::optional<std::string> request_timestamp;
        boost::optional<std::string> request_url;
        boost::optional<std::string> edition;
        boost::optional<std::string> about;
        boost::optional<std::string> pdgid;
        boost::optional<std::string> description;
        boost::optional<Summaries> summaries;

        public:
        boost::optional<int64_t> get_status_code() const { return status_code; }
        void set_status_code(boost::optional<int64_t> value) { this->status_code = value; }

        boost::optional<std::string> get_status_message() const { return status_message; }
        void set_status_message(boost::optional<std::string> value) { this->status_message = value; }

        boost::optional<std::string> get_request_timestamp() const { return request_timestamp; }
        void set_request_timestamp(boost::optional<std::string> value) { this->request_timestamp = value; }

        boost::optional<std::string> get_request_url() const { return request_url; }
        void set_request_url(boost::optional<std::string> value) { this->request_url = value; }

        boost::optional<std::string> get_edition() const { return edition; }
        void set_edition(boost::optional<std::string> value) { this->edition = value; }

        boost::optional<std::string> get_about() const { return about; }
        void set_about(boost::optional<std::string> value) { this->about = value; }

        boost::optional<std::string> get_pdgid() const { return pdgid; }
        void set_pdgid(boost::optional<std::string> value) { this->pdgid = value; }

        boost::optional<std::string> get_description() const { return description; }
        void set_description(boost::optional<std::string> value) { this->description = value; }

        boost::optional<Summaries> get_summaries() const { return summaries; }
        void set_summaries(boost::optional<Summaries> value) { this->summaries = value; }
    };
}

namespace PDGProvider {
    void from_json(const json & j, PdgValue & x);
    void to_json(json & j, const PdgValue & x);

    void from_json(const json & j, BranchingFraction & x);
    void to_json(json & j, const BranchingFraction & x);

    void from_json(const json & j, Property & x);
    void to_json(json & j, const Property & x);

    void from_json(const json & j, Summaries & x);
    void to_json(json & j, const Summaries & x);

    void from_json(const json & j, SummaryEdition & x);
    void to_json(json & j, const SummaryEdition & x);

    void from_json(const json & j, Type & x);
    void to_json(json & j, const Type & x);

    inline void from_json(const json & j, PdgValue& x) {
        x.set_value(j.at("value").get<double>());
        x.set_error_positive(get_stack_optional<double>(j, "error_positive"));
        x.set_error_negative(get_stack_optional<double>(j, "error_negative"));
        x.set_value_text(j.at("value_text").get<std::string>());
        x.set_type(j.at("type").get<Type>());
        x.set_scale_factor(get_stack_optional<double>(j, "scale_factor"));
        x.set_is_limit(get_stack_optional<bool>(j, "is_limit"));
        x.set_is_upper_limit(get_stack_optional<bool>(j, "is_upper_limit"));
        x.set_confidence_level(get_stack_optional<int64_t>(j, "confidence_level"));
        x.set_unit(get_stack_optional<std::string>(j, "unit"));
        x.set_comment(get_stack_optional<std::string>(j, "comment"));
    }

    inline void to_json(json & j, const PdgValue & x) {
        j = json::object();
        j["value"] = x.get_value();
        j["error_positive"] = x.get_error_positive();
        j["error_negative"] = x.get_error_negative();
        j["value_text"] = x.get_value_text();
        j["type"] = x.get_type();
        j["scale_factor"] = x.get_scale_factor();
        j["is_limit"] = x.get_is_limit();
        j["is_upper_limit"] = x.get_is_upper_limit();
        j["confidence_level"] = x.get_confidence_level();
        j["unit"] = x.get_unit();
        j["comment"] = x.get_comment();
    }

    inline void from_json(const json & j, BranchingFraction& x) {
        x.set_pdgid(get_stack_optional<std::string>(j, "pdgid"));
        x.set_description(get_stack_optional<std::string>(j, "description"));
        x.set_mode_number(get_stack_optional<int64_t>(j, "mode_number"));
        x.set_section(get_stack_optional<std::string>(j, "section"));
        x.set_pdg_values(get_stack_optional<std::vector<PdgValue>>(j, "pdg_values"));
    }

    inline void to_json(json & j, const BranchingFraction & x) {
        j = json::object();
        j["pdgid"] = x.get_pdgid();
        j["description"] = x.get_description();
        j["mode_number"] = x.get_mode_number();
        j["section"] = x.get_section();
        j["pdg_values"] = x.get_pdg_values();
    }

    inline void from_json(const json & j, Property& x) {
        x.set_pdgid(get_stack_optional<std::string>(j, "pdgid"));
        x.set_description(get_stack_optional<std::string>(j, "description"));
        x.set_pdg_values(get_stack_optional<std::vector<PdgValue>>(j, "pdg_values"));
    }

    inline void to_json(json & j, const Property & x) {
        j = json::object();
        j["pdgid"] = x.get_pdgid();
        j["description"] = x.get_description();
        j["pdg_values"] = x.get_pdg_values();
    }

    inline void from_json(const json & j, Summaries& x) {
        x.set_properties(get_stack_optional<std::vector<Property>>(j, "properties"));
        x.set_branching_fractions(get_stack_optional<std::vector<BranchingFraction>>(j, "branching_fractions"));
    }

    inline void to_json(json & j, const Summaries & x) {
        j = json::object();
        j["properties"] = x.get_properties();
        j["branching_fractions"] = x.get_branching_fractions();
    }

    inline void from_json(const json & j, SummaryEdition& x) {
        x.set_status_code(get_stack_optional<int64_t>(j, "status_code"));
        x.set_status_message(get_stack_optional<std::string>(j, "status_message"));
        x.set_request_timestamp(get_stack_optional<std::string>(j, "request_timestamp"));
        x.set_request_url(get_stack_optional<std::string>(j, "request_url"));
        x.set_edition(get_stack_optional<std::string>(j, "edition"));
        x.set_about(get_stack_optional<std::string>(j, "about"));
        x.set_pdgid(get_stack_optional<std::string>(j, "pdgid"));
        x.set_description(get_stack_optional<std::string>(j, "description"));
        x.set_summaries(get_stack_optional<Summaries>(j, "summaries"));
    }

    inline void to_json(json & j, const SummaryEdition & x) {
        j = json::object();
        j["status_code"] = x.get_status_code();
        j["status_message"] = x.get_status_message();
        j["request_timestamp"] = x.get_request_timestamp();
        j["request_url"] = x.get_request_url();
        j["edition"] = x.get_edition();
        j["about"] = x.get_about();
        j["pdgid"] = x.get_pdgid();
        j["description"] = x.get_description();
        j["summaries"] = x.get_summaries();
    }

    inline void from_json(const json & j, Type & x) {
        if (j == "BEST LIMIT") x = Type::BEST_LIMIT;
        else if (j == "OUR AVERAGE") x = Type::OUR_AVERAGE;
        else if (j == "OUR FIT") x = Type::OUR_FIT;
        else { throw std::runtime_error("Input JSON does not conform to schema!"); }
    }

    inline void to_json(json & j, const Type & x) {
        switch (x) {
            case Type::BEST_LIMIT: j = "BEST LIMIT"; break;
            case Type::OUR_AVERAGE: j = "OUR AVERAGE"; break;
            case Type::OUR_FIT: j = "OUR FIT"; break;
            default: throw std::runtime_error("Unexpected value in enumeration \"Type\": " + std::to_string(static_cast<int>(x)));
        }
    }
}
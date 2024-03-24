#include "AbstractSolver.h"

std::ostream& operator <<(std::ostream &stream, const std::vector<double> &vec){
    stream << "[";
    for (size_t i = 0 ;  i < vec.size(); i ++){
        stream << vec[i];
        if (i + 1 <vec.size()){
            stream << ", ";
        }
    }
    stream << "]";
    return stream;
}



void AbstractSolver::start_logging(size_t source, size_t target) {
    // All logging is done in JSON format
    std::stringstream start_info_json;
    start_info_json
        << "{\n"
        <<      "\t\"name\": \"" << get_solver_name() << "\",\n"
        <<      "\t\"eps\": " << this->eps << "\n"
        << "}";

    if (this->logger != nullptr) {
        LOG_START_SEARCH(*this->logger, source, target, start_info_json.str());
    }
}


void AbstractSolver::end_logging(SolutionSet &solutions, bool succ)
{
    // All logging is done in JSON format
    std::stringstream finish_info_json;
    finish_info_json
        << "{\n"
        << "\t\"solutions\": [";

    size_t solutions_count = 0;
    for (auto solution = solutions.begin(); solution != solutions.end(); ++solution) {
        if (solution != solutions.begin()) {
            finish_info_json << ",";
        }
        finish_info_json << "\n\t\t" << **solution;
        solutions_count++;
    }

    finish_info_json
        << "\n\t],\n"
        << "\t\"amount_of_solutions\": " << solutions_count << ",\n"
        << "\t\"status\": " << (succ ? "\"Success\"" : "\"Failed\"") << "\n"
        << "}" << std::endl;

    if (this->logger != nullptr) {
        LOG_FINISH_SEARCH(*(this->logger), finish_info_json.str());
    }
}

void AbstractSolver::log_optimal_paths(SolutionSet& solutions, std::string filename)
{
    std::cout << "Number of optimal paths = " << solutions.size() << std::endl;
    for (int i = 0; i < solutions.size(); i++)
    {
        std::cout << "====================\nPath " << (i + 1) << "\n====================\n";
        std::vector<Node> path;
        NodePtr currentNode = solutions[i];
        while (true)
        {
            path.push_back(*currentNode);
            if (currentNode->parent == nullptr)
            {
                break;
            }
            else
            {
                currentNode = currentNode->parent;
            }
        }
        for (int j = path.size() - 1; j >= 0; j--)
        {
            std::cout << "(" << (path.size() - j) << ") " << path[j].id << 
                " g(1) = " << path[j].g[0] << " g(2) = " << path[j].g[1] << std::endl;
        }
    }
    
    /*
    std::stringstream finish_info_json;
    finish_info_json
        << "{\n"
        << "\t\"solutions\": [";

    size_t solutions_count = 0;
    for (auto solution = solutions.begin(); solution != solutions.end(); ++solution) {
        if (solution != solutions.begin()) {
            finish_info_json << ",";
        }
        finish_info_json << "\n\t\t" << **solution;
        solutions_count++;
    }

    finish_info_json
        << "\n\t],\n"
        << "\t\"amount_of_solutions\": " << solutions_count << ",\n"
        << "\t\"number_of_expansions\": " << n_expansions << ",\n"
        << "\t\"number_of_generations\": " << n_generations << ",\n"
        << "\t\"status\": " << (succ ? "\"Success\"" : "\"Failed\"") << "\n"
        << "}" << std::endl;

    if (this->logger != nullptr) {
        LOG_FINISH_SEARCH(*(this->logger), finish_info_json.str());
    }
    */
}

void AbstractSolver::end_logging(SolutionSet & solutions, bool succ, int n_expansions, int n_generations)
{
    // All logging is done in JSON format
    std::stringstream finish_info_json;
    finish_info_json
        << "{\n"
        <<      "\t\"solutions\": [";

    size_t solutions_count = 0;
    for (auto solution = solutions.begin(); solution != solutions.end(); ++solution) {
        if (solution != solutions.begin()) {
            finish_info_json << ",";
        }
        finish_info_json << "\n\t\t" << **solution;
        solutions_count++;
    }

    finish_info_json
        <<      "\n\t],\n"
        <<      "\t\"amount_of_solutions\": " << solutions_count << ",\n"
        <<      "\t\"number_of_expansions\": " << n_expansions<< ",\n"
        <<      "\t\"number_of_generations\": " << n_generations << ",\n"
        <<      "\t\"status\": " << ( succ ? "\"Success\"": "\"Failed\"" )<< "\n"
        << "}" <<std::endl;

    if (this->logger != nullptr) {
        LOG_FINISH_SEARCH(*(this->logger), finish_info_json.str());
    }
}
